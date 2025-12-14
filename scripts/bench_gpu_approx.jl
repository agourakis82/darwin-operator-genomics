#!/usr/bin/env julia
"""
Benchmark GPU vs CPU performance for approximate symmetry computation.

Measures wall time and throughput for varying window counts and sizes.
"""

using BioSequences: LongDNA, randdnaseq
using Statistics
using Printf
using Dates

# Import GPU module
include(joinpath(@__DIR__, "..", "src", "gpu", "ApproximateSymmetryCUDA.jl"))
using .ApproximateSymmetryCUDA

function bench_dmin(windows::Vector, backend::String, n_reps::Int=3)::NamedTuple
    """Benchmark d_min computation with repetitions."""

    times = Float64[]

    for _ in 1:n_reps
        t_start = time()
        result = dmin_gpu(windows; backend=backend, include_rc=true)
        t_end = time()
        push!(times, (t_end - t_start) * 1000)  # Convert to ms
    end

    mean_time = mean(times)
    std_time = std(times)
    throughput = (length(windows) / (mean_time / 1000))  # windows per second

    return (
        mean_ms = mean_time,
        std_ms = std_time,
        min_ms = minimum(times),
        max_ms = maximum(times),
        throughput = throughput
    )
end

function main()
    println("\n=== GPU vs CPU Benchmark ===")
    println("Approximate Symmetry (d_min) Computation")
    println("Date: $(now())\n")

    # Test configurations
    window_lengths = [100, 250, 500, 1000]
    window_counts = [100, 500, 1000, 5000]

    println("Generating random test sequences...")
    println("Window Lengths: $window_lengths bp")
    println("Window Counts: $window_counts\n")

    # Create result table header
    println("=" ^ 100)
    @printf("%-8s %-10s %-12s %-12s %-12s %-15s\n",
            "L (bp)", "Count", "CPU Time(ms)", "CPU Tput(/s)", "CUDA Avail?", "Notes")
    println("=" ^ 100)

    cuda_available = has_cuda()

    for L in window_lengths
        for n_windows in window_counts
            # Generate random windows
            windows = [randdnaseq(L) for _ in 1:n_windows]

            # Bench CPU
            try
                cpu_result = bench_dmin(windows, "cpu", 2)

                # Bench CUDA if available
                cuda_result = nothing
                cuda_str = "No"
                if cuda_available
                    try
                        cuda_result = bench_dmin(windows, "cuda", 2)
                        cuda_str = "Yes"
                    catch e
                        cuda_str = "Error"
                    end
                end

                # Calculate speedup if both available
                speedup_str = ""
                if cuda_result !== nothing
                    speedup = cpu_result.mean_ms / cuda_result.mean_ms
                    speedup_str = @sprintf("%.2fx", speedup)
                end

                @printf("%-8d %-10d %-12.3f %-12.0f %-15s %s\n",
                        L,
                        n_windows,
                        cpu_result.mean_ms,
                        cpu_result.throughput,
                        cuda_str,
                        speedup_str)

                flush(stdout)
            catch e
                @printf("%-8d %-10d ERROR: %s\n", L, n_windows, string(e)[1:50])
            end
        end
        println("-" ^ 100)
    end

    println("\n=== Benchmark Summary ===")
    if cuda_available
        println("✓ CUDA.jl is available and tested")
        println("Note: Actual GPU speedup depends on:")
        println("  - GPU transfer overhead (PCIe bandwidth)")
        println("  - GPU kernel launch overhead")
        println("  - Problem size (larger windows benefit more)")
    else
        println("✗ CUDA.jl not available (CPU-only benchmark)")
        println("To enable GPU benchmarking, install CUDA.jl:")
        println("  julia> using Pkg; Pkg.add(\"CUDA\")")
    end

    println("\nPerformance Notes:")
    println("- CPU: Wall time for d_min computation across all windows")
    println("- Throughput: Windows per second")
    println("- Speedup: CPU time / GPU time (>1.0 = GPU faster)")
    println("- Typical GPU benefits emerge at 1000+ windows or very large L")

    println("\n=== Benchmark Complete ===")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
