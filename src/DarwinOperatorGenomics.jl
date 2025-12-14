module DarwinOperatorGenomics

using BioSequences
using Random
using StatsBase

export
    # Operators
    shift, reverse_seq, complement_seq, rev_comp,
    point_mutation, deletion, insertion, segment_inversion,
    # Canonicalization
    canonical_rep, orbit_size, is_fixed_under_reverse, is_fixed_under_rc,
    has_shift_symmetry,
    # Approximate symmetry
    min_dihedral_distance, hamming_distance, gc_shuffle,
    # Analysis
    compute_symmetry_stats, SymmetryStats,
    # Types
    GenomicOperator, OperatorChain, apply_chain,
    ShiftOp, ReverseOp, ComplementOp, ReverseComplementOp,
    MutationOp, DeletionOp, InsertionOp, InversionOp,
    operator_to_token, NUM_OPERATOR_TYPES, generate_random_chain

# =============================================================================
# Core Operators on DNA Sequences
# =============================================================================

"""
    shift(seq::LongDNA, k::Int) -> LongDNA

Cyclic shift: move first k bases to end. For circular genomes.
S^k(seq) = seq[k+1:end] * seq[1:k]
"""
function shift(seq::LongDNA, k::Int)
    n = length(seq)
    k = mod(k, n)
    k == 0 && return copy(seq)
    return seq[k+1:n] * seq[1:k]
end

"""
    reverse_seq(seq::LongDNA) -> LongDNA

Reverse operator R: reverses the sequence order.
R(seq) = seq[n:-1:1]
"""
function reverse_seq(seq::LongDNA)
    return BioSequences.reverse(seq)
end

"""
    complement_seq(seq::LongDNA) -> LongDNA

Complement operator K: A<->T, C<->G.
"""
function complement_seq(seq::LongDNA)
    return BioSequences.complement(seq)
end

"""
    rev_comp(seq::LongDNA) -> LongDNA

Reverse complement RC = R∘K = K∘R.
Aliased to avoid collision with BioSequences.reverse_complement.
"""
function rev_comp(seq::LongDNA)
    return BioSequences.reverse_complement(seq)
end

"""
    point_mutation(seq::LongDNA, pos::Int, newbase::DNA) -> LongDNA

M_{p->x}: Point substitution at position p to base x.
"""
function point_mutation(seq::LongDNA, pos::Int, newbase::DNA)
    result = copy(seq)
    result[pos] = newbase
    return result
end

"""
    deletion(seq::LongDNA, pos::Int) -> LongDNA

D_p: Delete base at position p. Non-invertible.
"""
function deletion(seq::LongDNA, pos::Int)
    n = length(seq)
    if pos == 1
        return seq[2:n]
    elseif pos == n
        return seq[1:n-1]
    else
        return seq[1:pos-1] * seq[pos+1:n]
    end
end

"""
    insertion(seq::LongDNA, pos::Int, base::DNA) -> LongDNA

I_{p,x}: Insert base x after position p. Non-invertible.
pos=0 inserts at beginning.
"""
function insertion(seq::LongDNA, pos::Int, base::DNA)
    n = length(seq)
    inserted = LongDNA{4}([base])
    if pos == 0
        return inserted * seq
    elseif pos >= n
        return seq * inserted
    else
        return seq[1:pos] * inserted * seq[pos+1:n]
    end
end

"""
    segment_inversion(seq::LongDNA, a::Int, b::Int) -> LongDNA

V_{a:b}: Invert segment from a to b (reverse complement the segment in place).
"""
function segment_inversion(seq::LongDNA, a::Int, b::Int)
    a, b = minmax(a, b)
    n = length(seq)
    inverted = BioSequences.reverse_complement(seq[a:b])
    if a == 1 && b == n
        return inverted
    elseif a == 1
        return inverted * seq[b+1:n]
    elseif b == n
        return seq[1:a-1] * inverted
    else
        return seq[1:a-1] * inverted * seq[b+1:n]
    end
end

# =============================================================================
# Canonicalization and Symmetry Analysis
# =============================================================================

"""
    canonical_rep(seq::LongDNA) -> LongDNA

Return lexicographically minimal representative under dihedral action
(all shifts and reverse of all shifts).
"""
function canonical_rep(seq::LongDNA)
    n = length(seq)
    n == 0 && return seq

    best = seq
    rev = reverse_seq(seq)

    for k in 0:n-1
        shifted = shift(seq, k)
        if shifted < best
            best = shifted
        end
        shifted_rev = shift(rev, k)
        if shifted_rev < best
            best = shifted_rev
        end
    end

    return best
end

"""
    orbit_size(seq::LongDNA) -> Int

Size of orbit under dihedral group action D_n.
Maximum is 2n (n shifts × 2 for reverse).
"""
function orbit_size(seq::LongDNA)
    n = length(seq)
    n == 0 && return 0

    seen = Set{LongDNA{4}}()
    rev = reverse_seq(seq)

    for k in 0:n-1
        push!(seen, shift(seq, k))
        push!(seen, shift(rev, k))
    end

    return length(seen)
end

"""
    is_fixed_under_reverse(seq::LongDNA) -> Bool

Check if seq is a palindrome (fixed under R).
"""
function is_fixed_under_reverse(seq::LongDNA)
    return seq == reverse_seq(seq)
end

"""
    is_fixed_under_rc(seq::LongDNA) -> Bool

Check if seq equals its reverse complement.
"""
function is_fixed_under_rc(seq::LongDNA)
    return seq == reverse_complement(seq)
end

"""
    has_shift_symmetry(seq::LongDNA) -> Tuple{Bool, Int}

Check if sequence has rotational symmetry.
Returns (has_symmetry, period) where period is smallest k>0 such that S^k(seq)=seq.
"""
function has_shift_symmetry(seq::LongDNA)
    n = length(seq)
    for k in 1:n-1
        if shift(seq, k) == seq
            return (true, k)
        end
    end
    return (false, n)
end

# =============================================================================
# Approximate Symmetry (d_min metric)
# =============================================================================

"""
    hamming_distance(a::LongDNA, b::LongDNA) -> Int

Count number of positions where sequences differ.
Sequences must have equal length.
"""
function hamming_distance(a::LongDNA, b::LongDNA)
    @assert length(a) == length(b) "Sequences must have equal length"
    d = 0
    for i in 1:length(a)
        if a[i] != b[i]
            d += 1
        end
    end
    return d
end

"""
    min_dihedral_distance(seq::LongDNA; include_rc::Bool=true) -> Int

Compute minimum Hamming distance to any non-identity dihedral transform.
d_min(w) = min over g ∈ {S^k, R∘S^k} \\ {id} of Hamming(w, g(w))

If include_rc=true, also includes RC∘S^k transforms.

Returns d_min as integer; divide by length(seq) for normalized distance.
"""
function min_dihedral_distance(seq::LongDNA; include_rc::Bool=true)
    n = length(seq)
    n == 0 && return 0

    d_min = n  # Initialize to max possible
    rev = reverse_seq(seq)

    # Check all shifts (excluding k=0 which is identity)
    for k in 1:n-1
        shifted = shift(seq, k)
        d = hamming_distance(seq, shifted)
        d_min = min(d_min, d)
    end

    # Check all reverse-shifts (including k=0)
    for k in 0:n-1
        shifted_rev = shift(rev, k)
        d = hamming_distance(seq, shifted_rev)
        d_min = min(d_min, d)
    end

    # Optionally check reverse-complement shifts
    if include_rc
        rc = rev_comp(seq)
        for k in 0:n-1
            shifted_rc = shift(rc, k)
            d = hamming_distance(seq, shifted_rc)
            d_min = min(d_min, d)
        end
    end

    return d_min
end

"""
    gc_shuffle(seq::LongDNA; rng=Random.GLOBAL_RNG) -> LongDNA

GC-preserving shuffle: randomly permute bases while preserving composition.
"""
function gc_shuffle(seq::LongDNA; rng=Random.GLOBAL_RNG)
    bases = collect(seq)
    shuffle!(rng, bases)
    return LongDNA{4}(bases)
end

# =============================================================================
# Operator Chains (for quaternion experiment)
# =============================================================================

abstract type GenomicOperator end

struct ShiftOp <: GenomicOperator
    k::Int
end

struct ReverseOp <: GenomicOperator end
struct ComplementOp <: GenomicOperator end
struct ReverseComplementOp <: GenomicOperator end

struct MutationOp <: GenomicOperator
    pos::Int
    newbase::DNA
end

struct DeletionOp <: GenomicOperator
    pos::Int
end

struct InsertionOp <: GenomicOperator
    pos::Int
    base::DNA
end

struct InversionOp <: GenomicOperator
    a::Int
    b::Int
end

const OperatorChain = Vector{GenomicOperator}

"""
    apply_op(seq::LongDNA, op::GenomicOperator) -> LongDNA

Apply a single operator to a sequence.
"""
function apply_op(seq::LongDNA, op::ShiftOp)
    shift(seq, op.k)
end

function apply_op(seq::LongDNA, op::ReverseOp)
    reverse_seq(seq)
end

function apply_op(seq::LongDNA, op::ComplementOp)
    complement_seq(seq)
end

function apply_op(seq::LongDNA, op::ReverseComplementOp)
    rev_comp(seq)
end

function apply_op(seq::LongDNA, op::MutationOp)
    point_mutation(seq, op.pos, op.newbase)
end

function apply_op(seq::LongDNA, op::DeletionOp)
    deletion(seq, op.pos)
end

function apply_op(seq::LongDNA, op::InsertionOp)
    insertion(seq, op.pos, op.base)
end

function apply_op(seq::LongDNA, op::InversionOp)
    segment_inversion(seq, op.a, op.b)
end

"""
    apply_chain(seq::LongDNA, chain) -> LongDNA

Apply a sequence of operators left-to-right.
Accepts any iterable of GenomicOperator subtypes.
"""
function apply_chain(seq::LongDNA, chain)
    result = seq
    for op in chain
        result = apply_op(result, op)
    end
    return result
end

"""
    operator_to_token(op::GenomicOperator) -> Int

Convert operator to integer token for ML models.
"""
function operator_to_token(op::GenomicOperator)
    if op isa ShiftOp
        return 1
    elseif op isa ReverseOp
        return 2
    elseif op isa ComplementOp
        return 3
    elseif op isa ReverseComplementOp
        return 4
    elseif op isa MutationOp
        return 5
    elseif op isa DeletionOp
        return 6
    elseif op isa InsertionOp
        return 7
    elseif op isa InversionOp
        return 8
    else
        return 0
    end
end

const NUM_OPERATOR_TYPES = 8

"""
    generate_random_chain(seq_len::Int, chain_len::Int; rng=Random.GLOBAL_RNG, include_indels=false) -> OperatorChain

Generate a random operator chain for simulation.
"""
function generate_random_chain(seq_len::Int, chain_len::Int;
                               rng=Random.GLOBAL_RNG, include_indels::Bool=false)
    bases = [DNA_A, DNA_C, DNA_G, DNA_T]
    chain = GenomicOperator[]
    current_len = seq_len

    for _ in 1:chain_len
        if include_indels
            op_type = rand(rng, 1:8)
        else
            op_type = rand(rng, 1:5)  # Exclude deletion, insertion, inversion
        end

        if op_type == 1  # Shift
            k = rand(rng, 1:current_len-1)
            push!(chain, ShiftOp(k))
        elseif op_type == 2  # Reverse
            push!(chain, ReverseOp())
        elseif op_type == 3  # Complement
            push!(chain, ComplementOp())
        elseif op_type == 4  # RC
            push!(chain, ReverseComplementOp())
        elseif op_type == 5  # Mutation
            pos = rand(rng, 1:current_len)
            base = rand(rng, bases)
            push!(chain, MutationOp(pos, base))
        elseif op_type == 6 && current_len > 10  # Deletion
            pos = rand(rng, 1:current_len)
            push!(chain, DeletionOp(pos))
            current_len -= 1
        elseif op_type == 7  # Insertion
            pos = rand(rng, 0:current_len)
            base = rand(rng, bases)
            push!(chain, InsertionOp(pos, base))
            current_len += 1
        elseif op_type == 8 && current_len > 10  # Inversion
            a = rand(rng, 1:current_len-5)
            b = rand(rng, a+1:min(a+100, current_len))
            push!(chain, InversionOp(a, b))
        else
            # Fallback to shift
            k = rand(rng, 1:max(1, current_len-1))
            push!(chain, ShiftOp(k))
        end
    end

    return chain
end

# =============================================================================
# Batch Analysis Functions
# =============================================================================

"""
    SymmetryStats

Results from symmetry analysis of a sequence.
"""
struct SymmetryStats
    length::Int
    orbit_size::Int
    max_orbit::Int  # 2n
    is_palindrome::Bool
    is_rc_fixed::Bool
    shift_period::Int
end

"""
    compute_symmetry_stats(seq::LongDNA) -> SymmetryStats
"""
function compute_symmetry_stats(seq::LongDNA)
    n = length(seq)
    orb = orbit_size(seq)
    has_sym, period = has_shift_symmetry(seq)

    return SymmetryStats(
        n,
        orb,
        2 * n,
        is_fixed_under_reverse(seq),
        is_fixed_under_rc(seq),
        period
    )
end

end # module
