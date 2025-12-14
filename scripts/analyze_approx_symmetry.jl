#!/usr/bin/env julia
"""
Analyze approximate dihedral self-similarity in bacterial genomes.

Computes d_min/L (minimum normalized distance to any dihedral transform)
across multiple window sizes, with GC-preserving shuffled baseline.
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
using Dates
using Random
using HypothesisTests

using DarwinOperatorGenomics

function parse_args()
    s = ArgParseSettings(
        description = "Analyze approximate dihedral symmetry",
        prog = "analyze_approx_symmetry.jl"
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
        "--windows", "-w"
            help = "Comma-separated window sizes"
            arg_type = String
            default = "100,250,500,1000"
        "--samples-per-size", "-n"
            help = "Windows to sample per genome per size"
            arg_type = Int
            default = 50
        "--seed", "-s"
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

"""
Sample windows and compute d_min/L for real and shuffled sequences.
"""
function analyze_windows(seq::LongDNA, window_sizes::Vector{Int}, n_samples::Int; rng=nothing)
    if rng === nothing
        rng = MersenneTwister(42)
    end

    results = DataFrame(
        window_size = Int[],
        d_min_normalized_real = Float64[],
        d_min_normalized_shuffled = Float64[]
    )

    n = length(seq)

    for L in window_sizes
        L >= n && continue

        max_start = n - L + 1
        positions = rand(rng, 1:max_start, min(n_samples, max_start))

        for pos in positions
            window = seq[pos:pos+L-1]

            # Real sequence d_min
            d_real = min_dihedral_distance(window; include_rc=true)

            # Shuffled baseline
            shuffled = gc_shuffle(window; rng=rng)
            d_shuffled = min_dihedral_distance(shuffled; include_rc=true)

            push!(results, (L, d_real / L, d_shuffled / L))
        end
    end

    return results
end

function main()
    args = parse_args()

    cache_dir = args["cache"]
    out_dir = args["out"]
    window_sizes = parse.(Int, split(args["windows"], ","))
    n_samples = args["samples-per-size"]
    seed = args["seed"]

    rng = MersenneTwister(seed)

    mkpath(joinpath(out_dir, "tables"))
    mkpath(joinpath(out_dir, "text"))

    println("\n=== Approximate Symmetry Analysis ===")
    println("Window sizes: $window_sizes")
    println("Samples per size per genome: $n_samples")

    manifest = load_manifest(cache_dir)
    println("Found $(length(manifest)) genomes")

    all_results = DataFrame(
        window_size = Int[],
        d_min_normalized_real = Float64[],
        d_min_normalized_shuffled = Float64[]
    )

    p = Progress(length(manifest); desc="Analyzing: ", showspeed=true)

    for entry in manifest
        filepath = entry["local_path"]
        !isfile(filepath) && (next!(p); continue)

        try
            sequences = read_fasta_gz(filepath)
            for (name, seq) in sequences
                length(seq) < maximum(window_sizes) * 2 && continue
                df = analyze_windows(seq, window_sizes, n_samples; rng=rng)
                append!(all_results, df)
            end
        catch e
            @warn "Error: $e"
        end
        next!(p)
    end

    # Save raw results
    csv_path = joinpath(out_dir, "tables", "approx_symmetry_stats.csv")
    CSV.write(csv_path, all_results)
    println("\nSaved: $csv_path")

    # Summary by window size
    summary_df = DataFrame(
        window_size = Int[],
        mean_d_min_real = Float64[],
        std_d_min_real = Float64[],
        mean_d_min_shuffled = Float64[],
        std_d_min_shuffled = Float64[],
        ks_pvalue = Float64[],
        effect_size = Float64[]
    )

    for L in window_sizes
        subset = filter(row -> row.window_size == L, all_results)
        nrow(subset) == 0 && continue

        real_vals = subset.d_min_normalized_real
        shuf_vals = subset.d_min_normalized_shuffled

        # KS test
        ks_result = ApproximateTwoSampleKSTest(real_vals, shuf_vals)
        pval = pvalue(ks_result)

        # Effect size (Cohen's d)
        pooled_std = sqrt((var(real_vals) + var(shuf_vals)) / 2)
        cohens_d = pooled_std > 0 ? (mean(real_vals) - mean(shuf_vals)) / pooled_std : 0.0

        push!(summary_df, (
            L,
            mean(real_vals),
            std(real_vals),
            mean(shuf_vals),
            std(shuf_vals),
            pval,
            cohens_d
        ))
    end

    summary_path = joinpath(out_dir, "tables", "approx_symmetry_summary.csv")
    CSV.write(summary_path, summary_df)
    println("Saved: $summary_path")

    # Text report
    text = """
# Results: Approximate Dihedral Self-Similarity

## Method

For each window of length L, we compute:
- **d_min/L**: Minimum normalized Hamming distance to any dihedral transform g(w), excluding identity
- **Baseline**: Same metric on GC-preserving shuffled sequences

## Summary by Window Size

| L (bp) | d_min/L (Real) | d_min/L (Shuffled) | KS p-value | Effect Size |
|--------|----------------|--------------------|-----------:|------------:|
"""

    for row in eachrow(summary_df)
        text *= "| $(row.window_size) | $(round(row.mean_d_min_real, digits=4)) ± $(round(row.std_d_min_real, digits=4)) | $(round(row.mean_d_min_shuffled, digits=4)) ± $(round(row.std_d_min_shuffled, digits=4)) | $(round(row.ks_pvalue, sigdigits=3)) | $(round(row.effect_size, digits=3)) |\n"
    end

    text *= """

## Interpretation

- **d_min/L < 0.5**: Strong approximate self-similarity under some dihedral transform
- **Real ≈ Shuffled**: No excess structure beyond base composition
- **Real < Shuffled**: Genomic sequences have more dihedral self-similarity than random
- **Effect size > 0.2**: Small practical difference; > 0.5: medium; > 0.8: large

---
*Generated: $(now())*
"""

    text_path = joinpath(out_dir, "text", "results_approx_symmetry.md")
    write(text_path, text)
    println("Saved: $text_path")

    println("\n=== Analysis Complete ===")
end

main()
