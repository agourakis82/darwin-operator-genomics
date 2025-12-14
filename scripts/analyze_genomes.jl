#!/usr/bin/env julia
"""
Analyze bacterial genomes for dihedral symmetry statistics.

Reads downloaded FASTA files and computes:
- Orbit sizes under dihedral group action
- Fixed points under R and RC
- Summary statistics and distributions
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

# Load our module
using DarwinOperatorGenomics

function parse_args()
    s = ArgParseSettings(
        description = "Analyze bacterial genomes for symmetry statistics",
        prog = "analyze_genomes.jl"
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
            help = "Window size for local symmetry analysis"
            arg_type = Int
            default = 1000
        "--sample-windows", "-n"
            help = "Number of windows to sample per genome"
            arg_type = Int
            default = 100
        "--seed", "-s"
            help = "Random seed"
            arg_type = Int
            default = 42
    end

    return ArgParse.parse_args(s)
end

"""
Load manifest and return list of genome files.
"""
function load_manifest(cache_dir::String)
    manifest_path = joinpath(cache_dir, "manifest.jsonl")

    if !isfile(manifest_path)
        error("Manifest not found: $manifest_path. Run fetch_ncbi_complete_bacteria.jl first.")
    end

    entries = Dict[]
    for line in eachline(manifest_path)
        isempty(strip(line)) && continue
        push!(entries, JSON3.read(line, Dict))
    end

    return entries
end

"""
Read all sequences from a gzipped FASTA file.
Returns vector of (name, sequence) tuples.
"""
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
Analyze a single replicon for symmetry statistics.
"""
function analyze_replicon(seq::LongDNA, window_size::Int, n_windows::Int; rng=nothing)
    n = length(seq)

    # For very short sequences, analyze the whole thing
    if n <= window_size * 2
        stats = compute_symmetry_stats(seq)
        return DataFrame(
            length = [n],
            orbit_size = [stats.orbit_size],
            max_orbit = [stats.max_orbit],
            orbit_ratio = [stats.orbit_size / stats.max_orbit],
            is_palindrome = [stats.is_palindrome],
            is_rc_fixed = [stats.is_rc_fixed],
            shift_period = [stats.shift_period],
            has_shift_symmetry = [stats.shift_period < n]
        )
    end

    # Sample random windows
    results = DataFrame(
        length = Int[],
        orbit_size = Int[],
        max_orbit = Int[],
        orbit_ratio = Float64[],
        is_palindrome = Bool[],
        is_rc_fixed = Bool[],
        shift_period = Int[],
        has_shift_symmetry = Bool[]
    )

    # Deterministic window positions
    if rng === nothing
        rng = Random.MersenneTwister(42)
    end

    max_start = n - window_size + 1
    positions = sort(unique(rand(rng, 1:max_start, n_windows)))

    for pos in positions
        window = seq[pos:pos+window_size-1]
        stats = compute_symmetry_stats(window)

        push!(results, (
            window_size,
            stats.orbit_size,
            stats.max_orbit,
            stats.orbit_size / stats.max_orbit,
            stats.is_palindrome,
            stats.is_rc_fixed,
            stats.shift_period,
            stats.shift_period < window_size
        ))
    end

    return results
end

"""
Classify replicon type from FASTA header.
"""
function classify_replicon(name::String)
    name_lower = lowercase(name)
    if occursin("plasmid", name_lower)
        return "plasmid"
    elseif occursin("chromosome", name_lower)
        return "chromosome"
    else
        return "unknown"
    end
end

function main()
    args = parse_args()

    cache_dir = args["cache"]
    out_dir = args["out"]
    window_size = args["window"]
    n_windows = args["sample-windows"]
    seed = args["seed"]

    rng = MersenneTwister(seed)

    mkpath(out_dir)
    mkpath(joinpath(out_dir, "tables"))
    mkpath(joinpath(out_dir, "text"))

    # Load manifest
    println("\n=== Loading Manifest ===")
    manifest = load_manifest(cache_dir)
    println("Found $(length(manifest)) genomes in manifest")

    # Collect results
    all_replicon_stats = DataFrame(
        accession = String[],
        organism = String[],
        replicon_name = String[],
        replicon_type = String[],
        replicon_length = Int[],
        mean_orbit_ratio = Float64[],
        palindrome_fraction = Float64[],
        rc_fixed_fraction = Float64[],
        shift_sym_fraction = Float64[]
    )

    all_window_stats = DataFrame(
        accession = String[],
        organism = String[],
        replicon_type = String[],
        window_length = Int[],
        orbit_size = Int[],
        max_orbit = Int[],
        orbit_ratio = Float64[],
        is_palindrome = Bool[],
        is_rc_fixed = Bool[]
    )

    println("\n=== Analyzing Genomes ===")
    p = Progress(length(manifest); desc="Analyzing: ", showspeed=true)

    for entry in manifest
        accession = entry["accession"]
        organism = entry["organism"]
        filepath = entry["local_path"]

        if !isfile(filepath)
            @warn "File not found: $filepath"
            next!(p)
            continue
        end

        # Read sequences
        try
            sequences = read_fasta_gz(filepath)

            for (name, seq) in sequences
                replicon_type = classify_replicon(name)
                n = length(seq)

                # Skip very short sequences
                n < 100 && continue

                # Analyze windows
                window_df = analyze_replicon(seq, window_size, n_windows; rng=rng)

                # Aggregate stats for this replicon
                if nrow(window_df) > 0
                    push!(all_replicon_stats, (
                        accession,
                        organism,
                        name,
                        replicon_type,
                        n,
                        mean(window_df.orbit_ratio),
                        mean(window_df.is_palindrome),
                        mean(window_df.is_rc_fixed),
                        mean(window_df.has_shift_symmetry)
                    ))

                    # Add window-level stats
                    for row in eachrow(window_df)
                        push!(all_window_stats, (
                            accession,
                            organism,
                            replicon_type,
                            row.length,
                            row.orbit_size,
                            row.max_orbit,
                            row.orbit_ratio,
                            row.is_palindrome,
                            row.is_rc_fixed
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
    println("\n=== Saving Results ===")

    # Table 1: Dataset summary
    table1 = DataFrame(
        metric = String[],
        value = Any[]
    )

    push!(table1, ("N assemblies", length(manifest)))
    push!(table1, ("N replicons analyzed", nrow(all_replicon_stats)))

    n_chr = count(==("chromosome"), all_replicon_stats.replicon_type)
    n_pla = count(==("plasmid"), all_replicon_stats.replicon_type)
    n_unk = count(==("unknown"), all_replicon_stats.replicon_type)

    push!(table1, ("N chromosomes", n_chr))
    push!(table1, ("N plasmids", n_pla))
    push!(table1, ("N unknown", n_unk))

    if nrow(all_replicon_stats) > 0
        push!(table1, ("Mean replicon length (bp)", round(Int, mean(all_replicon_stats.replicon_length))))
        push!(table1, ("Median replicon length (bp)", round(Int, median(all_replicon_stats.replicon_length))))
        push!(table1, ("Min replicon length (bp)", minimum(all_replicon_stats.replicon_length)))
        push!(table1, ("Max replicon length (bp)", maximum(all_replicon_stats.replicon_length)))
    end

    table1_path = joinpath(out_dir, "tables", "table1_dataset_summary.csv")
    CSV.write(table1_path, table1)
    println("Saved: $table1_path")

    # Save detailed stats
    replicon_stats_path = joinpath(out_dir, "tables", "replicon_stats.csv")
    CSV.write(replicon_stats_path, all_replicon_stats)
    println("Saved: $replicon_stats_path")

    window_stats_path = joinpath(out_dir, "tables", "window_stats.csv")
    CSV.write(window_stats_path, all_window_stats)
    println("Saved: $window_stats_path")

    # Generate text summary
    println("\n=== Generating Text Summary ===")

    text_summary = """
# Results: Dihedral Symmetry Analysis

## Dataset Summary

- **N assemblies**: $(length(manifest))
- **N replicons**: $(nrow(all_replicon_stats))
- **Chromosomes**: $n_chr
- **Plasmids**: $n_pla

## Symmetry Statistics (window size = $window_size bp)

"""

    if nrow(all_window_stats) > 0
        mean_orbit_ratio = mean(all_window_stats.orbit_ratio)
        std_orbit_ratio = std(all_window_stats.orbit_ratio)
        pal_frac = mean(all_window_stats.is_palindrome)
        rc_frac = mean(all_window_stats.is_rc_fixed)

        text_summary *= """
### Overall Statistics

- **Mean orbit ratio** (orbit_size / 2n): $(round(mean_orbit_ratio, digits=4)) ± $(round(std_orbit_ratio, digits=4))
- **Palindrome fraction**: $(round(pal_frac * 100, digits=2))%
- **RC-fixed fraction**: $(round(rc_frac * 100, digits=2))%

### Interpretation

An orbit ratio near 1.0 indicates the sequence has no dihedral symmetry (all rotations and reflections produce distinct sequences). Values < 1.0 indicate symmetry.

- Orbit ratio ≈ 1.0: No detectable symmetry at this window size
- Palindrome fraction: Proportion of windows that read the same forwards and backwards
- RC-fixed fraction: Proportion of windows equal to their reverse complement

"""

        # Compare chromosomes vs plasmids if both present
        if n_chr > 0 && n_pla > 0
            chr_stats = filter(row -> row.replicon_type == "chromosome", all_window_stats)
            pla_stats = filter(row -> row.replicon_type == "plasmid", all_window_stats)

            if nrow(chr_stats) > 0 && nrow(pla_stats) > 0
                text_summary *= """
### Chromosome vs Plasmid Comparison

| Metric | Chromosomes | Plasmids |
|--------|-------------|----------|
| Mean orbit ratio | $(round(mean(chr_stats.orbit_ratio), digits=4)) | $(round(mean(pla_stats.orbit_ratio), digits=4)) |
| Palindrome % | $(round(mean(chr_stats.is_palindrome)*100, digits=2)) | $(round(mean(pla_stats.is_palindrome)*100, digits=2)) |
| RC-fixed % | $(round(mean(chr_stats.is_rc_fixed)*100, digits=2)) | $(round(mean(pla_stats.is_rc_fixed)*100, digits=2)) |

"""
            end
        end
    end

    text_summary *= """
---
*Generated: $(now())*
"""

    text_path = joinpath(out_dir, "text", "results_r1_r2.md")
    write(text_path, text_summary)
    println("Saved: $text_path")

    println("\n=== Analysis Complete ===")
end

main()
