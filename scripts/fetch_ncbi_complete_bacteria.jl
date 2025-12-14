#!/usr/bin/env julia
"""
Fetch complete bacterial genomes from NCBI RefSeq.

Downloads assembly_summary.txt, filters for complete genomes,
and downloads FASTA files with checksums.
"""

using ArgParse
using HTTP
using CSV
using DataFrames
using CodecZlib
using SHA
using JSON3
using Random
using ProgressMeter
using Dates

function parse_args()
    s = ArgParseSettings(
        description = "Download NCBI complete bacterial genomes",
        prog = "fetch_ncbi_complete_bacteria.jl"
    )

    @add_arg_table! s begin
        "--max", "-m"
            help = "Maximum number of assemblies to download"
            arg_type = Int
            default = 200
        "--out", "-o"
            help = "Output directory for cache"
            arg_type = String
            default = "data/cache/"
        "--seed", "-s"
            help = "Random seed for reproducible sampling"
            arg_type = Int
            default = 42
        "--resume"
            help = "Resume from existing manifest"
            action = :store_true
        "--dry-run"
            help = "List assemblies without downloading"
            action = :store_true
    end

    return ArgParse.parse_args(s)
end

const ASSEMBLY_SUMMARY_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt"

"""
Download and parse NCBI assembly summary.
"""
function fetch_assembly_summary(cache_dir::String)
    summary_path = joinpath(cache_dir, "assembly_summary.txt")

    if !isfile(summary_path)
        println("Downloading assembly_summary.txt...")
        response = HTTP.get(ASSEMBLY_SUMMARY_URL; readtimeout=120)
        write(summary_path, response.body)
        println("  Saved to $summary_path")
    else
        println("Using cached assembly_summary.txt")
    end

    # Parse TSV (skip comment lines starting with ##)
    lines = readlines(summary_path)
    header_idx = findfirst(l -> startswith(l, "#assembly_accession"), lines)

    if header_idx === nothing
        error("Could not find header line in assembly_summary.txt")
    end

    # Extract header (remove leading "#")
    header_line = lines[header_idx][2:end]
    headers = split(header_line, '\t')

    # Read data lines
    data_lines = filter(l -> !startswith(l, "#"), lines[header_idx+1:end])

    # Build DataFrame
    df = DataFrame([Symbol(h) => String[] for h in headers])

    for line in data_lines
        fields = split(line, '\t')
        if length(fields) >= length(headers)
            push!(df, fields[1:length(headers)])
        end
    end

    return df
end

"""
Filter for complete genomes with latest version.
"""
function filter_complete_genomes(df::DataFrame)
    # Column names from NCBI assembly_summary
    # assembly_level is column 12 (index 12), version_status is column 11

    # Find column indices
    level_col = :assembly_level
    status_col = :version_status

    filtered = filter(row ->
        row[level_col] == "Complete Genome" &&
        row[status_col] == "latest",
        df
    )

    println("Found $(nrow(filtered)) complete genomes with 'latest' status")
    return filtered
end

"""
Compute SHA256 checksum of a file.
"""
function compute_sha256(filepath::String)
    return bytes2hex(open(sha256, filepath))
end

"""
Download a single genome FASTA.
"""
function download_genome(ftp_path::String, accession::String, cache_dir::String)
    # Convert FTP to HTTPS
    https_base = replace(ftp_path, "ftp://" => "https://")

    # Build filename
    basename_part = split(ftp_path, '/')[end]
    fasta_filename = "$(basename_part)_genomic.fna.gz"
    fasta_url = "$https_base/$(fasta_filename)"

    local_path = joinpath(cache_dir, "genomes", fasta_filename)

    if isfile(local_path)
        return (path=local_path, url=fasta_url, cached=true)
    end

    mkpath(dirname(local_path))

    try
        response = HTTP.get(fasta_url; readtimeout=300)
        write(local_path, response.body)
        return (path=local_path, url=fasta_url, cached=false)
    catch e
        @warn "Failed to download $accession: $e"
        return nothing
    end
end

"""
Write manifest entry to JSONL file.
"""
function write_manifest_entry(manifest_path::String, entry::Dict)
    open(manifest_path, "a") do io
        JSON3.write(io, entry)
        println(io)
    end
end

"""
Load existing manifest entries.
"""
function load_manifest(manifest_path::String)
    existing = Dict{String, Dict}()

    if isfile(manifest_path)
        for line in eachline(manifest_path)
            isempty(strip(line)) && continue
            entry = JSON3.read(line, Dict)
            existing[entry["accession"]] = entry
        end
    end

    return existing
end

function main()
    args = parse_args()

    cache_dir = args["out"]
    max_assemblies = args["max"]
    seed = args["seed"]
    resume = args["resume"]
    dry_run = args["dry-run"]

    mkpath(cache_dir)
    mkpath(joinpath(cache_dir, "genomes"))

    manifest_path = joinpath(cache_dir, "manifest.jsonl")

    # Fetch and filter assembly summary
    println("\n=== Fetching NCBI Assembly Summary ===")
    df = fetch_assembly_summary(cache_dir)
    println("Total assemblies in summary: $(nrow(df))")

    complete = filter_complete_genomes(df)

    # Random sampling if needed
    if nrow(complete) > max_assemblies
        Random.seed!(seed)
        indices = randperm(nrow(complete))[1:max_assemblies]
        complete = complete[indices, :]
        println("Sampled $max_assemblies assemblies (seed=$seed)")
    end

    if dry_run
        println("\n=== Dry Run: Would download $(nrow(complete)) assemblies ===")
        for row in eachrow(complete)
            println("  $(row[:assembly_accession]): $(row[:organism_name])")
        end
        return
    end

    # Load existing manifest for resume
    existing = resume ? load_manifest(manifest_path) : Dict{String, Dict}()
    println("Existing manifest entries: $(length(existing))")

    # Download genomes
    println("\n=== Downloading Genomes ===")
    downloaded = 0
    skipped = 0
    failed = 0

    p = Progress(nrow(complete); desc="Downloading: ", showspeed=true)

    for row in eachrow(complete)
        accession = row[:assembly_accession]
        ftp_path = row[:ftp_path]
        organism = row[:organism_name]

        if haskey(existing, accession)
            skipped += 1
            next!(p)
            continue
        end

        result = download_genome(ftp_path, accession, cache_dir)

        if result === nothing
            failed += 1
            next!(p)
            continue
        end

        # Compute checksum
        sha = compute_sha256(result.path)
        filesize = stat(result.path).size

        # Write manifest entry
        entry = Dict(
            "accession" => accession,
            "organism" => organism,
            "url" => result.url,
            "local_path" => result.path,
            "timestamp" => string(now()),
            "size_bytes" => filesize,
            "sha256" => sha,
            "cached" => result.cached
        )

        write_manifest_entry(manifest_path, entry)
        downloaded += 1
        next!(p)
    end

    println("\n=== Summary ===")
    println("Downloaded: $downloaded")
    println("Skipped (existing): $skipped")
    println("Failed: $failed")
    println("Manifest: $manifest_path")
end

main()
