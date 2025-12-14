#!/usr/bin/env julia
"""
Analyze inverted repeats enrichment in bacterial genomes.

Detects IRs and compares to Markov-shuffled baseline.
"""

using ArgParse
using JSON3
using FASTX
using BioSequences
using CodecZlib
using DataFrames
using CSV
using Statistics
using Random
using ProgressMeter

include(joinpath(@__DIR__, "..", "src", "InvertedRepeats.jl"))
using .InvertedRepeats

function parse_args()
    s = ArgParseSettings(
        description = "Analyze inverted repeats enrichment",
        prog = "analyze_inverted_repeats.jl"
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
        "--min-stem"
            help = "Minimum stem length"
            arg_type = Int
            default = 8
        "--max-loop"
            help = "Maximum loop length"
            arg_type = Int
            default = 20
        "--min-loop"
            help = "Minimum loop length"
            arg_type = Int
            default = 3
        "--n-samples"
            help = "Number of baseline samples"
            arg_type = Int
            default = 5
        "--max-genomes"
            help = "Max genomes to analyze (0 for all)"
            arg_type = Int
            default = 0
        "--seed"
            help = "Random seed"
            arg_type = Int
            default = 42
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
    min_stem = args["min-stem"]
    max_loop = args["max-loop"]
    min_loop = args["min-loop"]
    n_samples = args["n-samples"]
    max_genomes = args["max-genomes"]
    seed = args["seed"]

    rng = MersenneTwister(seed)

    mkpath(joinpath(out_dir, "tables"))

    println("\n=== Inverted Repeats Enrichment Analysis ===")
    println("Stem: ≥$min_stem, Loop: $min_loop-$max_loop, Samples: $n_samples")

    manifest = load_manifest(cache_dir)
    println("Found $(length(manifest)) genomes")

    # Limit if specified
    if max_genomes > 0 && length(manifest) > max_genomes
        manifest = manifest[1:max_genomes]
        println("Limited to $max_genomes genomes")
    end

    # Results storage
    enrichment_df = DataFrame(
        replicon_id = String[],
        replicon_length = Int[],
        observed_ir_count = Int[],
        baseline_mean = Float64[],
        baseline_std = Float64[],
        enrichment_ratio = Float64[],
        z_score = Float64[]
    )

    # Detailed IR list (sample)
    ir_details_df = DataFrame(
        replicon_id = String[],
        position = Int[],
        stem_length = Int[],
        loop_length = Int[],
        stem_sequence = String[]
    )

    detail_limit = 1000  # Limit detailed IR storage

    p = Progress(length(manifest); desc="Analyzing: ", showspeed=true)

    for entry in manifest
        filepath = entry["local_path"]
        !isfile(filepath) && (next!(p); continue)

        try
            sequences = read_fasta_gz(filepath)
            for (name, seq) in sequences
                n = length(seq)
                n < 50000 && continue  # Skip small sequences

                # Enrichment test
                result = ir_enrichment_test(
                    seq, name;
                    min_stem=min_stem, max_loop=max_loop, min_loop=min_loop,
                    n_samples=n_samples, shuffle_order=1, rng=rng
                )

                push!(enrichment_df, (
                    name,
                    n,
                    result.observed_count,
                    result.baseline_mean,
                    result.baseline_std,
                    result.enrichment_ratio,
                    result.z_score
                ))

                # Store detailed IRs for first few genomes
                if nrow(ir_details_df) < detail_limit
                    irs = find_inverted_repeats(seq; min_stem=min_stem, max_loop=max_loop, min_loop=min_loop)
                    for ir in irs
                        nrow(ir_details_df) >= detail_limit && break
                        push!(ir_details_df, (
                            name,
                            ir.pos,
                            ir.stem_len,
                            ir.loop_len,
                            ir.stem_seq
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
    enrichment_path = joinpath(out_dir, "tables", "ir_enrichment_summary.csv")
    CSV.write(enrichment_path, enrichment_df)
    println("\nSaved: $enrichment_path")

    details_path = joinpath(out_dir, "tables", "inverted_repeats.csv")
    CSV.write(details_path, ir_details_df)
    println("Saved: $details_path")

    # Summary statistics
    println("\n=== Summary ===")
    println("Total replicons analyzed: $(nrow(enrichment_df))")

    if nrow(enrichment_df) > 0
        println("\nObserved IR counts:")
        println("  median: $(median(enrichment_df.observed_ir_count))")
        println("  mean: $(round(mean(enrichment_df.observed_ir_count), digits=1))")
        println("  max: $(maximum(enrichment_df.observed_ir_count))")

        println("\nEnrichment ratio (observed/baseline):")
        finite_ratios = filter(isfinite, enrichment_df.enrichment_ratio)
        if !isempty(finite_ratios)
            println("  median: $(round(median(finite_ratios), digits=2))")
            println("  mean: $(round(mean(finite_ratios), digits=2))")
        end

        println("\nZ-scores:")
        finite_z = filter(isfinite, enrichment_df.z_score)
        if !isempty(finite_z)
            println("  median: $(round(median(finite_z), digits=2))")
            println("  mean: $(round(mean(finite_z), digits=2))")
            enriched = count(z -> z > 2.0, finite_z)
            println("  enriched (z > 2): $enriched / $(length(finite_z)) ($(round(100*enriched/length(finite_z), digits=1))%)")
        end
    end

    println("\n=== Analysis Complete ===")
end

main()
