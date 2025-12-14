"""
K-mer Inversion Symmetry Module

Quantifies generalized Chargaff / inversion symmetry: the degree to which
k-mer counts equal their reverse-complement counts on a single strand.

Reference: Albrecht-Buehler (2006), Forsdyke & Mortimer (2000), BMC Genomics 2016.
"""
module KmerInversionSymmetry

using BioSequences
using Statistics

export
    count_kmers, kmer_inversion_score, compute_inversion_symmetry_profile,
    find_symmetry_limit, InversionSymmetryResult

"""
    count_kmers(seq::LongDNA, k::Int) -> Dict{LongDNA{4}, Int}

Count all k-mers in sequence. Skips k-mers containing ambiguous bases.
"""
function count_kmers(seq::LongDNA, k::Int)
    counts = Dict{LongDNA{4}, Int}()
    n = length(seq)
    k > n && return counts

    for i in 1:(n - k + 1)
        kmer = seq[i:i+k-1]
        # Skip if contains N or other ambiguous
        has_ambig = false
        for j in 1:k
            b = kmer[j]
            if !(b == DNA_A || b == DNA_C || b == DNA_G || b == DNA_T)
                has_ambig = true
                break
            end
        end
        has_ambig && continue
        counts[kmer] = get(counts, kmer, 0) + 1
    end
    return counts
end

"""
    kmer_inversion_score(counts::Dict{LongDNA{4}, Int}, k::Int) -> Float64

Compute inversion symmetry score X_k:
    X_k = mean over unique k-mer pairs {w, RC(w)} of |N(w) - N(RC(w))| / (N(w) + N(RC(w)) + ε)

Returns value in [0, 1]: 0 = perfect symmetry, 1 = maximal asymmetry.
"""
function kmer_inversion_score(counts::Dict{LongDNA{4}, Int}, k::Int)
    seen = Set{LongDNA{4}}()
    scores = Float64[]
    eps = 1e-10

    for (kmer, n_w) in counts
        kmer in seen && continue
        rc_kmer = BioSequences.reverse_complement(kmer)
        push!(seen, kmer)
        push!(seen, rc_kmer)

        n_rc = get(counts, rc_kmer, 0)

        # Asymmetry score for this pair
        score = abs(n_w - n_rc) / (n_w + n_rc + eps)
        push!(scores, score)
    end

    isempty(scores) && return 0.0
    return mean(scores)
end

"""
    InversionSymmetryResult

Result of inversion symmetry analysis for a single replicon.
"""
struct InversionSymmetryResult
    replicon_id::String
    replicon_length::Int
    k_values::Vector{Int}
    X_k::Vector{Float64}           # Inversion score per k
    K_L_tau01::Int                  # K-limit at τ=0.1
    K_L_tau005::Int                 # K-limit at τ=0.05
end

"""
    compute_inversion_symmetry_profile(seq::LongDNA, replicon_id::String; K_max::Int=10) -> InversionSymmetryResult

Compute inversion symmetry scores X_k for k=1..K_max.
"""
function compute_inversion_symmetry_profile(seq::LongDNA, replicon_id::String; K_max::Int=10)
    n = length(seq)
    k_values = Int[]
    X_k = Float64[]

    for k in 1:min(K_max, n)
        counts = count_kmers(seq, k)
        score = kmer_inversion_score(counts, k)
        push!(k_values, k)
        push!(X_k, score)
    end

    # Find K_L for different tolerances
    K_L_01 = find_symmetry_limit(X_k, 0.1)
    K_L_005 = find_symmetry_limit(X_k, 0.05)

    return InversionSymmetryResult(
        replicon_id,
        n,
        k_values,
        X_k,
        K_L_01,
        K_L_005
    )
end

"""
    find_symmetry_limit(X_k::Vector{Float64}, tau::Float64) -> Int

Find K_L: the largest k such that X_k ≤ τ.
Returns 0 if no k satisfies the condition.
"""
function find_symmetry_limit(X_k::Vector{Float64}, tau::Float64)
    K_L = 0
    for (k, x) in enumerate(X_k)
        if x <= tau
            K_L = k
        end
    end
    return K_L
end

"""
    compute_inversion_symmetry_by_half(seq::LongDNA, replicon_id::String, ori_pos::Int; K_max::Int=10)

Compute inversion symmetry separately for each replichore (half of circular genome).
ori_pos is 1-indexed position of replication origin.
Returns tuple of (leading_half_result, lagging_half_result).
"""
function compute_inversion_symmetry_by_half(seq::LongDNA, replicon_id::String, ori_pos::Int; K_max::Int=10)
    n = length(seq)
    half_len = div(n, 2)

    # Split at ori: leading strand goes ori -> ter, lagging goes ter -> ori
    # For circular: rotate so ori is at position 1
    rotated = seq[ori_pos:end] * seq[1:ori_pos-1]

    leading = rotated[1:half_len]
    lagging = rotated[half_len+1:end]

    result_leading = compute_inversion_symmetry_profile(leading, "$(replicon_id)_leading"; K_max=K_max)
    result_lagging = compute_inversion_symmetry_profile(lagging, "$(replicon_id)_lagging"; K_max=K_max)

    return (result_leading, result_lagging)
end

end # module
