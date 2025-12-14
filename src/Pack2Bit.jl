"""
    Pack2Bit

GPU-friendly 2-bit DNA encoding for CUDA acceleration of symmetry computations.

Encodes DNA bases as 2-bit values:
- A (Adenine)   → 0
- C (Cytosine)  → 1
- G (Guanine)   → 2
- T (Thymine)   → 3

Multiple bases are packed into UInt64 chunks for efficient XOR + popcount operations.
"""
module Pack2Bit

using BioSequences: LongDNA, DNA_A, DNA_C, DNA_G, DNA_T

export pack_window_2bit, unpack_window_2bit, chunks_needed, pack_multiple_windows, base_to_bits, bits_to_base

"""
    base_to_bits(base) -> UInt8

Convert a DNA base to 2-bit representation (0-3).
- A → 0
- C → 1
- G → 2
- T → 3
"""
function base_to_bits(base)::UInt8
    if base == DNA_A
        return 0x00
    elseif base == DNA_C
        return 0x01
    elseif base == DNA_G
        return 0x02
    elseif base == DNA_T
        return 0x03
    else
        error("Unknown DNA base: $base")
    end
end

"""
    bits_to_base(bits::UInt8)

Convert 2-bit value back to DNA base.
"""
function bits_to_base(bits::UInt8)
    bits = bits & 0x03  # Mask to 2 bits
    if bits == 0x00
        return DNA_A
    elseif bits == 0x01
        return DNA_C
    elseif bits == 0x02
        return DNA_G
    else  # bits == 0x03
        return DNA_T
    end
end

"""
    chunks_needed(length::Int) -> Int

Calculate number of UInt64 chunks needed to store a sequence of given length.
Each UInt64 holds 32 bases (32 * 2 bits = 64 bits).
"""
function chunks_needed(length::Int)::Int
    return div(length + 31, 32)  # Ceiling division
end

"""
    pack_window_2bit(seq::LongDNA) -> Vector{UInt64}

Pack a DNA sequence into 2-bit chunks stored in UInt64 array.

Each UInt64 stores 32 bases with 2 bits per base.
Bases are stored little-endian within each chunk:
- Bits [1:2]    → base 1
- Bits [3:4]    → base 2
- ...
- Bits [63:64]  → base 32

Returns a vector of UInt64 chunks.
"""
function pack_window_2bit(seq::LongDNA)::Vector{UInt64}
    L = length(seq)
    n_chunks = chunks_needed(L)
    packed = zeros(UInt64, n_chunks)

    for i = 1:L
        base = seq[i]
        bits = base_to_bits(base)
        chunk_idx = div(i - 1, 32) + 1  # Which chunk (1-indexed)
        bit_offset = (mod(i - 1, 32) * 2)  # Bit position within chunk

        packed[chunk_idx] |= (UInt64(bits) << bit_offset)
    end

    return packed
end

"""
    unpack_window_2bit(packed::Vector{UInt64}, length::Int)

Unpack a 2-bit packed window back to DNA sequence.

Args:
- packed: Vector of UInt64 chunks
- length: Original sequence length (to trim padding)
"""
function unpack_window_2bit(packed::Vector{UInt64}, length::Int)
    bases = []

    for i = 1:length
        chunk_idx = div(i - 1, 32) + 1
        bit_offset = (mod(i - 1, 32) * 2)

        bits = UInt8((packed[chunk_idx] >> bit_offset) & 0x03)
        push!(bases, bits_to_base(bits))
    end

    return LongDNA{4}(bases)
end

"""
    pack_multiple_windows(windows::Vector) -> (Matrix{UInt64}, Int)

Pack multiple windows into a matrix for GPU transfer.

Returns:
- packed: Matrix of shape (n_chunks, n_windows) containing packed bases
- L: Length of windows (assumed uniform)
"""
function pack_multiple_windows(windows::Vector)
    if isempty(windows)
        error("Cannot pack empty window list")
    end

    L = length(windows[1])
    for w in windows
        if length(w) != L
            error("All windows must have same length")
        end
    end

    n_windows = length(windows)
    n_chunks = chunks_needed(L)

    packed = zeros(UInt64, n_chunks, n_windows)

    for (w_idx, window) in enumerate(windows)
        packed[:, w_idx] = pack_window_2bit(window)
    end

    return packed, L
end

end  # module Pack2Bit
