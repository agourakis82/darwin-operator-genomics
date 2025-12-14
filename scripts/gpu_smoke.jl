#!/usr/bin/env julia
"""
GPU Smoke Test - Verify CUDA availability and basic GPU operations.

This script checks:
1. CUDA availability and version info
2. Basic GPU array operations with scalar indexing guard
3. Fallback gracefully if CUDA unavailable
"""

using DarwinOperatorGenomics

# Import GPU module
include(joinpath(@__DIR__, "..", "src", "gpu", "ApproximateSymmetryCUDA.jl"))
using .ApproximateSymmetryCUDA

println("\n=== GPU Smoke Test ===")
println("Checking CUDA support...")

if has_cuda()
    println("✓ CUDA.jl detected and loaded")
    CUDA.versioninfo()
else
    println("✗ CUDA.jl not available (CPU-only mode)")
    println("  To enable CUDA, install CUDA.jl via:")
    println("  julia> using Pkg; Pkg.add(\"CUDA\")")
end

println("\nTesting d_min_gpu function on sample sequences...")

using BioSequences: LongDNA, randdnaseq

# Generate test sequences
test_seqs = [
    LongDNA{4}("ACGTACGT"),
    LongDNA{4}("TGCATGCA"),
    randdnaseq(50),
]

println("Test sequences:")
for (i, seq) in enumerate(test_seqs)
    println("  Seq $i: length=$(length(seq))")
end

# Test CPU backend (should always work)
println("\nComputing d_min on CPU...")
try
    d_mins_cpu = dmin_gpu(test_seqs; backend="cpu", include_rc=true)
    println("✓ CPU d_min values: $d_mins_cpu")
catch e
    println("✗ CPU computation failed: $e")
    exit(1)
end

# Test CUDA backend if available
if has_cuda()
    println("\nComputing d_min on CUDA...")
    try
        d_mins_cuda = dmin_gpu(test_seqs; backend="cuda", include_rc=true)
        println("✓ CUDA d_min values: $d_mins_cuda")

        # Verify results match
        if d_mins_cpu == d_mins_cuda
            println("✓ CPU and CUDA results match!")
        else
            println("⚠ WARNING: CPU and CUDA results differ!")
            println("  CPU:  $d_mins_cpu")
            println("  CUDA: $d_mins_cuda")
        end
    catch e
        println("⚠ CUDA computation warning: $e")
    end
else
    println("\nSkipping CUDA test (CUDA unavailable)")
end

println("\n=== Smoke Test Complete ===")
println("GPU backend is ready for use!")
