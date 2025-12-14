#!/usr/bin/env julia
"""
Analyze GC skew and estimate ori/ter positions for bacterial genomes.

Computes cumulative GC skew to estimate replication geometry,
then stratifies symmetry metrics by replichore.
"""

using ArgParse
using JSON3
using FASTX
using BioSequences
using CodecZlib
using DataFrames
using CSV
using Statistics
using ProgressMeter

include(joinpath(@__DIR__, "..", "src", "GCSkew.jl"))
include(joinpath(@__DIR__, "..", "src", "KmerInversionSymmetry.jl"))
using .GCSkew
using .KmerInversionSymmetry

using DarwinOperatorGenomics

function parse_args()
    s = ArgParseSettings(
        description = "Analyze GC skew and replichore symmetry",
        prog = "analyze_gc_skew.jl"
    )

    @add_arg_table! s begin
        "--cache", "-c"
            help = "Cache directory with downloaded genomes"
            arg_type = String
            default = "data/cache/"
        "--out", "-o"
            help = "Output directory for results"
            arg_type = String
            default = "results/"
        "--window", "-w"
            help = "Window size for GC skew"
            arg_type = Int
            default = 1000
        "--step", "-s"
            help = "Step size for GC skew"
            arg_type = Int
            default = 500
        "--kmax"
            help = "Max k for inversion symmetry"
            arg_type = Int
            default = 10
    end

    return ArgParse.parse_args(s)
end

function load_manifest(cache_dir::String)
    manifest_path = joinpath(cache_dir, "manifest.jsonl")
    entries = Dict[]
    if isfile(manifest_path)
        for line in eachline(manifest_path)
            isempty(strip(line)) && continue
            push!(entries, JSON3.read(line, Dict))
        end
    end
    return entries
end

function read_fasta_gz(filepath::String)
    sequences = Tuple{String, LongDNA{4}}[]
    open(GzipDecompressorStream, filepath) do io
        reader = FASTA.Reader(io)
        for record in reader
            name = FASTA.identifier(record)
            seq = LongDNA{4}(FASTA.sequence(record))
            push!(sequences, (name, seq))
        end
    end
    return sequences
end

function main()
    args = parse_args()

    cache_dir = args["cache"]
    out_dir = args["out"]
    window_size = args["window"]
    step_size = args["step"]
    K_max = args["kmax"]

    mkpath(joinpath(out_dir, "tables"))

    println("\n=== GC Skew & Replichore Analysis ===")
    println("Window: $window_size, Step: $step_size")

    manifest = load_manifest(cache_dir)
    println("Found $(length(manifest)) genomes")

    # Results storage
    ori_ter_df = DataFrame(
        replicon_id = String[],
        replicon_length = Int[],
        ori_pos = Int[],
        ter_pos = Int[],
        skew_amplitude = Float64[],
        confidence = String[],
        leading_gc = Float64[],
        lagging_gc = Float64[],
        gc_diff = Float64[]
    )

    replichore_df = DataFrame(
        replicon_id = String[],
        replichore = String[],
        length = Int[],
        X_1 = Float64[],
        X_5 = Float64[],
        X_10 = Float64[],
        K_L_tau01 = Int[],
        K_L_tau005 = Int[]
    )

    p = Progress(length(manifest); desc="Analyzing: ", showspeed=true)

    for entry in manifest
        filepath = entry["local_path"]
        !isfile(filepath) && (next!(p); continue)

        try
            sequences = read_fasta_gz(filepath)
            for (name, seq) in sequences
                n = length(seq)
                n < 50000 && continue  # Need substantial length for skew

                # Estimate ori/ter
                estimate = estimate_ori_ter(seq; window_size=window_size, step=step_size)

                # Replichore GC content
                gc_stats = compute_replichore_gc_content(seq, estimate.ori_pos)

                push!(ori_ter_df, (
                    name,
                    n,
                    estimate.ori_pos,
                    estimate.ter_pos,
                    estimate.skew_amplitude,
                    string(estimate.confidence),
                    gc_stats.leading_gc,
                    gc_stats.lagging_gc,
                    gc_stats.gc_diff
                ))

                # Only do replichore stratification for high/medium confidence
                if estimate.confidence in [:high, :medium]
                    leading, lagging = split_replichores(seq, estimate.ori_pos)

                    for (rep_seq, rep_name) in [(leading, "leading"), (lagging, "lagging")]
                        length(rep_seq) < 10000 && continue

                        result = compute_inversion_symmetry_profile(
                            rep_seq,
                            "$(name)_$(rep_name)";
                            K_max=K_max
                        )

                        X_1 = length(result.X_k) >= 1 ? result.X_k[1] : NaN
                        X_5 = length(result.X_k) >= 5 ? result.X_k[5] : NaN
                        X_10 = length(result.X_k) >= 10 ? result.X_k[10] : NaN

                        push!(replichore_df, (
                            name,
                            rep_name,
                            length(rep_seq),
                            X_1, X_5, X_10,
                            result.K_L_tau01,
                            result.K_L_tau005
                        ))
                    end
                end
            end
        catch e
            @warn "Error processing $filepath: $e"
        end
        next!(p)
    end

    # Save results
    ori_ter_path = joinpath(out_dir, "tables", "gc_skew_ori_ter.csv")
    CSV.write(ori_ter_path, ori_ter_df)
    println("\nSaved: $ori_ter_path")

    replichore_path = joinpath(out_dir, "tables", "replichore_symmetry.csv")
    CSV.write(replichore_path, replichore_df)
    println("Saved: $replichore_path")

    # Summary statistics
    println("\n=== Summary ===")
    println("Total replicons: $(nrow(ori_ter_df))")

    if nrow(ori_ter_df) > 0
        conf_counts = combine(groupby(ori_ter_df, :confidence), nrow => :count)
        println("\nConfidence distribution:")
        for row in eachrow(conf_counts)
            println("  $(row.confidence): $(row.count)")
        end

        println("\nSkew amplitude:")
        println("  median: $(round(median(ori_ter_df.skew_amplitude), digits=4))")
        println("  mean: $(round(mean(ori_ter_df.skew_amplitude), digits=4))")
    end

    if nrow(replichore_df) > 0
        println("\nReplichore comparison (high/medium confidence only):")

        leading = filter(r -> r.replichore == "leading", replichore_df)
        lagging = filter(r -> r.replichore == "lagging", replichore_df)

        if nrow(leading) > 0 && nrow(lagging) > 0
            println("  Leading X_5 median: $(round(median(skipmissing(leading.X_5)), digits=4))")
            println("  Lagging X_5 median: $(round(median(skipmissing(lagging.X_5)), digits=4))")
        end
    end

    println("\n=== Analysis Complete ===")
end

main()
