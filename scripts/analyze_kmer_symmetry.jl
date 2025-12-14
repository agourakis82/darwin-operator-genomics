#!/usr/bin/env julia
"""
Analyze k-mer inversion symmetry across bacterial genomes.

Computes generalized Chargaff symmetry: X_k scores for k=1..K_max,
and finds K_L (symmetry limit) at different tolerances.
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

include(joinpath(@__DIR__, "..", "src", "KmerInversionSymmetry.jl"))
using .KmerInversionSymmetry

function parse_args()
    s = ArgParseSettings(
        description = "Analyze k-mer inversion symmetry",
        prog = "analyze_kmer_symmetry.jl"
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
        "--kmax", "-k"
            help = "Maximum k-mer size"
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
    K_max = args["kmax"]

    mkpath(joinpath(out_dir, "tables"))

    println("\n=== K-mer Inversion Symmetry Analysis ===")
    println("K_max: $K_max")

    manifest = load_manifest(cache_dir)
    println("Found $(length(manifest)) genomes")

    # Results storage
    all_results = DataFrame(
        replicon_id = String[],
        replicon_length = Int[],
        k = Int[],
        X_k = Float64[]
    )

    summary_results = DataFrame(
        replicon_id = String[],
        replicon_length = Int[],
        K_L_tau01 = Int[],
        K_L_tau005 = Int[],
        X_1 = Float64[],
        X_5 = Float64[],
        X_10 = Float64[]
    )

    p = Progress(length(manifest); desc="Analyzing: ", showspeed=true)

    for entry in manifest
        filepath = entry["local_path"]
        !isfile(filepath) && (next!(p); continue)

        try
            sequences = read_fasta_gz(filepath)
            for (name, seq) in sequences
                length(seq) < 10000 && continue  # Skip very short

                result = compute_inversion_symmetry_profile(seq, name; K_max=K_max)

                # Per-k results
                for (k, xk) in zip(result.k_values, result.X_k)
                    push!(all_results, (name, result.replicon_length, k, xk))
                end

                # Summary
                X_1 = length(result.X_k) >= 1 ? result.X_k[1] : NaN
                X_5 = length(result.X_k) >= 5 ? result.X_k[5] : NaN
                X_10 = length(result.X_k) >= 10 ? result.X_k[10] : NaN

                push!(summary_results, (
                    name,
                    result.replicon_length,
                    result.K_L_tau01,
                    result.K_L_tau005,
                    X_1, X_5, X_10
                ))
            end
        catch e
            @warn "Error processing $filepath: $e"
        end
        next!(p)
    end

    # Save raw results
    csv_path = joinpath(out_dir, "tables", "kmer_inversion_symmetry.csv")
    CSV.write(csv_path, all_results)
    println("\nSaved: $csv_path")

    # Save summary
    summary_path = joinpath(out_dir, "tables", "kmer_inversion_symmetry_summary.csv")
    CSV.write(summary_path, summary_results)
    println("Saved: $summary_path")

    # Compute aggregates by k
    agg_df = DataFrame(
        k = Int[],
        median_X_k = Float64[],
        q25_X_k = Float64[],
        q75_X_k = Float64[],
        mean_X_k = Float64[],
        std_X_k = Float64[],
        n_replicons = Int[]
    )

    for k in 1:K_max
        subset = filter(row -> row.k == k, all_results)
        nrow(subset) == 0 && continue

        vals = subset.X_k
        push!(agg_df, (
            k,
            median(vals),
            quantile(vals, 0.25),
            quantile(vals, 0.75),
            mean(vals),
            std(vals),
            nrow(subset)
        ))
    end

    agg_path = joinpath(out_dir, "tables", "kmer_inversion_symmetry_by_k.csv")
    CSV.write(agg_path, agg_df)
    println("Saved: $agg_path")

    # Print summary
    println("\n=== Summary ===")
    println("Total replicons analyzed: $(nrow(summary_results))")
    if nrow(agg_df) > 0
        println("\nMedian X_k by k:")
        for row in eachrow(agg_df)
            println("  k=$(row.k): median=$(round(row.median_X_k, digits=4)), IQR=[$(round(row.q25_X_k, digits=4)), $(round(row.q75_X_k, digits=4))]")
        end
    end

    # K_L distribution
    if nrow(summary_results) > 0
        println("\nK_L distribution (τ=0.1):")
        println("  median: $(median(summary_results.K_L_tau01))")
        println("  mean: $(round(mean(summary_results.K_L_tau01), digits=1))")

        println("\nK_L distribution (τ=0.05):")
        println("  median: $(median(summary_results.K_L_tau005))")
        println("  mean: $(round(mean(summary_results.K_L_tau005), digits=1))")
    end

    println("\n=== Analysis Complete ===")
end

main()
