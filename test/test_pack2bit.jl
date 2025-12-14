using Test
using BioSequences: LongDNA, DNA_A, DNA_C, DNA_G, DNA_T, randdnaseq
include(joinpath(@__DIR__, "..", "src", "Pack2Bit.jl"))
using .Pack2Bit

@testset "Pack2Bit" begin
    @testset "base_to_bits / bits_to_base round-trip" begin
        for base in [DNA_A, DNA_C, DNA_G, DNA_T]
            bits = base_to_bits(base)
            recovered = bits_to_base(bits)
            @test recovered == base
        end
    end

    @testset "chunks_needed calculation" begin
        @test chunks_needed(1) == 1
        @test chunks_needed(32) == 1
        @test chunks_needed(33) == 2
        @test chunks_needed(64) == 2
        @test chunks_needed(65) == 3
    end

    @testset "pack_window_2bit / unpack_window_2bit round-trip" begin
        # Test short sequence
        seq_short = LongDNA{4}("ACGT")
        packed = pack_window_2bit(seq_short)
        unpacked = unpack_window_2bit(packed, length(seq_short))
        @test unpacked == seq_short

        # Test longer sequence with all bases repeated
        seq_long = LongDNA{4}("ACGTACGTACGTACGT")  # 16 bases
        packed = pack_window_2bit(seq_long)
        unpacked = unpack_window_2bit(packed, length(seq_long))
        @test unpacked == seq_long

        # Test sequence longer than 32 bases (requires multiple chunks)
        bases = "ACGT" * "ACGT"^8  # 32 bases: one full chunk
        seq_full = LongDNA{4}(bases)
        packed = pack_window_2bit(seq_full)
        unpacked = unpack_window_2bit(packed, length(seq_full))
        @test unpacked == seq_full

        # Test with more than one chunk
        bases_multi = "ACGT"^20  # 80 bases: 3 chunks
        seq_multi = LongDNA{4}(bases_multi)
        packed = pack_window_2bit(seq_multi)
        unpacked = unpack_window_2bit(packed, length(seq_multi))
        @test unpacked == seq_multi
    end

    @testset "pack_multiple_windows" begin
        windows = [
            LongDNA{4}("ACGTACGTAA"),
            LongDNA{4}("TGCATGCAAA"),
            LongDNA{4}("AAAAAAAAAA")
        ]
        packed, L = pack_multiple_windows(windows)
        @test size(packed, 2) == 3
        @test L == 10

        # Verify round-trip
        for i = 1:length(windows)
            unpacked = unpack_window_2bit(packed[:, i], L)
            @test unpacked == windows[i]
        end
    end

    @testset "random sequence packing" begin
        using Random
        Random.seed!(42)

        for len in [50, 100, 128, 256, 1000]
            random_seq = randdnaseq(len)
            packed = pack_window_2bit(random_seq)
            unpacked = unpack_window_2bit(packed, len)
            @test unpacked == random_seq
        end
    end
end
