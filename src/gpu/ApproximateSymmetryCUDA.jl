"""
    ApproximateSymmetryCUDA

Optional GPU acceleration for approximate dihedral symmetry computation.
Falls back to CPU if CUDA is not available.

This module provides GPU-accelerated computation of d_min (minimum Hamming distance
to dihedral transforms) for packed DNA windows.

To use: include this file and call `dmin_gpu()` with backend="cuda" or "cpu".
"""
module ApproximateSymmetryCUDA

using BioSequences: LongDNA, reverse, reverse_complement

export dmin_gpu, has_cuda, dmin_cpu, dmin_cuda_impl, hamming_distance

# Try to load CUDA.jl if available
const HAS_CUDA = try
    using CUDA
    true
catch
    false
end

if HAS_CUDA
    export CUDA
end

"""
    has_cuda() -> Bool

Check if CUDA support is available.
"""
function has_cuda()::Bool
    return HAS_CUDA
end

"""
    dmin_cpu(windows::Vector, include_rc::Bool) -> Vector{Int}

CPU implementation: compute d_min for each window.
Falls back when CUDA unavailable or explicitly requested.

Args:
- windows: Vector of LongDNA sequences (all same length)
- include_rc: Whether to test reverse-complement transforms

Returns:
- Vector of Int: d_min distance for each window
"""
function dmin_cpu(windows::Vector, include_rc::Bool)::Vector{Int}
    # Import CPU d_min function from main module
    # This avoids circular dependencies by computing d_min directly here

    d_mins = Int[]

    for window in windows
        L = length(window)
        d_min = L  # Max possible distance

        # S^k transforms (cyclic shifts)
        for k in 1:L-1
            shifted = window[k+1:end] * window[1:k]  # Julia string-like indexing
            d = hamming_distance(window, shifted)
            d_min = min(d_min, d)
        end

        # R∘S^k transforms (reverse then shift)
        rev = reverse(window)
        for k in 0:L-1
            shifted_rev = if k == 0
                rev
            else
                rev[k+1:end] * rev[1:k]
            end
            d = hamming_distance(window, shifted_rev)
            d_min = min(d_min, d)
        end

        # RC∘S^k transforms (reverse complement then shift)
        if include_rc
            rc = reverse_complement(window)
            for k in 0:L-1
                shifted_rc = if k == 0
                    rc
                else
                    rc[k+1:end] * rc[1:k]
                end
                d = hamming_distance(window, shifted_rc)
                d_min = min(d_min, d)
            end
        end

        push!(d_mins, d_min)
    end

    return d_mins
end

"""
    hamming_distance(a::LongDNA, b::LongDNA) -> Int

Compute Hamming distance (number of mismatches) between two DNA sequences.
"""
function hamming_distance(a::LongDNA, b::LongDNA)::Int
    @assert length(a) == length(b)
    d = 0
    for i in 1:length(a)
        if a[i] != b[i]
            d += 1
        end
    end
    return d
end

"""
    dmin_gpu(windows::Vector; backend::String="cpu", include_rc::Bool=true) -> Vector{Int}

Compute d_min for each window using specified backend.

Args:
- windows: Vector of DNA sequences (all same length)
- backend: "cpu" (default) or "cuda"
- include_rc: Whether to test reverse-complement transforms (default: true)

Returns:
- Vector of d_min values (one per window)

Note: If CUDA unavailable, falls back to CPU regardless of backend setting.
"""
function dmin_gpu(windows::Vector; backend::String="cpu", include_rc::Bool=true)::Vector{Int}
    if backend == "cuda" && HAS_CUDA
        return dmin_cuda_impl(windows, include_rc)
    else
        # CPU fallback (default or CUDA unavailable)
        return dmin_cpu(windows, include_rc)
    end
end

"""
    dmin_cuda_impl(windows::Vector{LongDNA}, include_rc::Bool) -> Vector{Int}

CUDA kernel implementation (requires CUDA.jl).
Placeholder for future GPU acceleration.
Currently falls back to CPU.

TODO (Phase 3): Implement packed 2-bit encoding and CUDA kernels via @cuda.
"""
function dmin_cuda_impl(windows::Vector{LongDNA}, include_rc::Bool)::Vector{Int}
    # Placeholder: Future implementation will:
    # 1. Pack windows into 2-bit representation (Pack2Bit module)
    # 2. Transfer to GPU via CuArray
    # 3. Run @cuda kernels to compute Hamming distances
    # 4. Return results

    # For now, fall back to CPU
    @warn "CUDA.jl available but GPU kernels not yet implemented. Using CPU."
    return dmin_cpu(windows, include_rc)
end

end  # module ApproximateSymmetryCUDA
