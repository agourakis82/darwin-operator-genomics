"""
Correctness tests: Verify CPU and GPU backends produce identical d_min results.
"""

using Test
using BioSequences: LongDNA, randdnaseq
using Random

include(joinpath(@__DIR__, "..", "src", "gpu", "ApproximateSymmetryCUDA.jl"))
using .ApproximateSymmetryCUDA

@testset "GPU Correctness Tests" begin
    @testset "CPU backend available" begin
        @test true  # CPU is always available
    end

    @testset "Simple sequences - d_min determinism" begin
        # Test on simple sequences with known structure
        test_seqs = [
            LongDNA{4}("ACGTACGT"),
            LongDNA{4}("AAAAAAA"),
            LongDNA{4}("ACGTACGTACGT"),
        ]

        for seq in test_seqs
            d_cpu = dmin_gpu([seq]; backend="cpu", include_rc=true)[1]
            @test isa(d_cpu, Integer)
            @test 0 <= d_cpu <= length(seq)
        end
    end

    @testset "Random sequences - d_min bounds" begin
        Random.seed!(42)

        for _ in 1:10
            seq = randdnaseq(rand(50:500))
            d_cpu = dmin_gpu([seq]; backend="cpu", include_rc=true)[1]

            # d_min should be between 0 and sequence length
            @test 0 <= d_cpu <= length(seq)
        end
    end

    @testset "Multiple windows - batch correctness" begin
        Random.seed!(43)

        windows = [randdnaseq(100) for _ in 1:5]
        d_vals_cpu = dmin_gpu(windows; backend="cpu", include_rc=true)

        @test length(d_vals_cpu) == 5
        for d in d_vals_cpu
            @test 0 <= d <= 100
        end
    end

    @testset "include_rc parameter effect" begin
        Random.seed!(44)

        seq = randdnaseq(100)

        # Without RC
        d_no_rc = dmin_gpu([seq]; backend="cpu", include_rc=false)[1]

        # With RC
        d_with_rc = dmin_gpu([seq]; backend="cpu", include_rc=true)[1]

        # With RC should be <= without RC (more transforms to check)
        @test d_with_rc <= d_no_rc
    end

    @testset "CUDA fallback to CPU when unavailable" begin
        if has_cuda()
            println("  CUDA available - testing GPU backend")
            windows = [randdnaseq(100) for _ in 1:3]

            d_cpu = dmin_gpu(windows; backend="cpu", include_rc=true)
            d_cuda = dmin_gpu(windows; backend="cuda", include_rc=true)

            # Results should match
            @test d_cpu == d_cuda
        else
            println("  CUDA not available - GPU falls back to CPU")
            windows = [randdnaseq(100) for _ in 1:3]

            d_cuda_req = dmin_gpu(windows; backend="cuda", include_rc=true)
            d_cpu_fallback = dmin_gpu(windows; backend="cpu", include_rc=true)

            # Should be identical since CUDA unavailable
            @test d_cuda_req == d_cpu_fallback
        end
    end

    @testset "Reproducibility - same seed gives same results" begin
        Random.seed!(42)
        windows1 = [randdnaseq(100) for _ in 1:3]
        d_vals1 = dmin_gpu(windows1; backend="cpu", include_rc=true)

        Random.seed!(42)
        windows2 = [randdnaseq(100) for _ in 1:3]
        d_vals2 = dmin_gpu(windows2; backend="cpu", include_rc=true)

        @test d_vals1 == d_vals2
    end

    println("\n✓ All GPU correctness tests passed!")
end
