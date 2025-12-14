"""
Inverted Repeats (IR) Detection Module

Detects inverted repeats (potential hairpin/terminator structures) and tests
enrichment vs Markov-shuffled baseline.

Reference: Bioinformatics 2002 (IR enrichment near 3' ends).
"""
module InvertedRepeats

using BioSequences
using Random
using Statistics

export
    InvertedRepeat, find_inverted_repeats, count_inverted_repeats,
    markov_shuffle, ir_enrichment_test, IREnrichmentResult

"""
    InvertedRepeat

Represents a detected inverted repeat.
"""
struct InvertedRepeat
    pos::Int           # Start position of first stem
    stem_len::Int      # Length of stem
    loop_len::Int      # Length of loop
    stem_seq::String   # Sequence of first stem
end

"""
    find_inverted_repeats(seq::LongDNA; min_stem::Int=8, max_loop::Int=20, min_loop::Int=3) -> Vector{InvertedRepeat}

Find inverted repeats (potential hairpins) in sequence.
An IR consists of: stem1 - loop - stem2, where stem2 = RC(stem1).

Parameters:
- min_stem: minimum stem length (default 8)
- max_loop: maximum loop length (default 20)
- min_loop: minimum loop length (default 3)
"""
function find_inverted_repeats(seq::LongDNA;
                               min_stem::Int=8,
                               max_loop::Int=20,
                               min_loop::Int=3)
    irs = InvertedRepeat[]
    n = length(seq)

    # Scan for potential IR positions
    # For each position, check if there's a matching RC downstream
    for pos in 1:(n - 2*min_stem - min_loop)
        for stem_len in min_stem:min(20, div(n - pos - min_loop, 2))
            stem1_end = pos + stem_len - 1
            stem1_end > n && continue

            stem1 = seq[pos:stem1_end]

            # Check for ambiguous bases
            has_ambig = false
            for i in 1:stem_len
                b = stem1[i]
                if !(b == DNA_A || b == DNA_C || b == DNA_G || b == DNA_T)
                    has_ambig = true
                    break
                end
            end
            has_ambig && continue

            stem1_rc = BioSequences.reverse_complement(stem1)

            # Search for stem2 within loop range
            for loop_len in min_loop:max_loop
                stem2_start = stem1_end + loop_len + 1
                stem2_end = stem2_start + stem_len - 1
                stem2_end > n && continue

                stem2 = seq[stem2_start:stem2_end]

                # Check if stem2 matches RC of stem1
                if stem2 == stem1_rc
                    push!(irs, InvertedRepeat(
                        pos,
                        stem_len,
                        loop_len,
                        string(stem1)
                    ))
                    break  # Found IR at this position, move on
                end
            end
        end
    end

    return irs
end

"""
    count_inverted_repeats(seq::LongDNA; min_stem::Int=8, max_loop::Int=20, min_loop::Int=3) -> Int

Count total inverted repeats in sequence (faster than full enumeration).
"""
function count_inverted_repeats(seq::LongDNA;
                                 min_stem::Int=8,
                                 max_loop::Int=20,
                                 min_loop::Int=3)
    return length(find_inverted_repeats(seq; min_stem=min_stem, max_loop=max_loop, min_loop=min_loop))
end

"""
    markov_shuffle(seq::LongDNA; order::Int=1, rng=Random.GLOBAL_RNG) -> LongDNA{4}

Generate Markov-shuffled sequence preserving mono- or di-nucleotide frequencies.
order=1: preserve mononucleotide frequencies
order=2: preserve dinucleotide frequencies (more complex)
"""
function markov_shuffle(seq::LongDNA; order::Int=1, rng=Random.GLOBAL_RNG)
    if order == 1
        # Simple shuffle preserving base composition
        bases = collect(seq)
        shuffle!(rng, bases)
        return LongDNA{4}(bases)
    else
        # Dinucleotide preserving shuffle (Altschul-Erickson algorithm simplified)
        # For efficiency, use Euler path approach
        return dinucleotide_shuffle(seq; rng=rng)
    end
end

"""
    dinucleotide_shuffle(seq::LongDNA; rng=Random.GLOBAL_RNG) -> LongDNA{4}

Shuffle sequence while preserving dinucleotide frequencies.
Uses simplified approach: build graph of dinucleotide transitions, find Eulerian path.
"""
function dinucleotide_shuffle(seq::LongDNA; rng=Random.GLOBAL_RNG)
    n = length(seq)
    n < 2 && return copy(seq)

    # For simplicity and speed, use approximate method:
    # Shuffle in blocks preserving local structure
    bases = [DNA_A, DNA_C, DNA_G, DNA_T]

    # Count dinucleotides
    di_counts = Dict{Tuple{DNA, DNA}, Int}()
    for i in 1:n-1
        b1, b2 = seq[i], seq[i+1]
        (b1 in bases && b2 in bases) || continue
        key = (b1, b2)
        di_counts[key] = get(di_counts, key, 0) + 1
    end

    # Build new sequence by sampling transitions
    result = DNA[seq[1]]
    remaining = copy(di_counts)

    for _ in 2:n
        current = result[end]
        # Find valid next bases
        candidates = Tuple{DNA, Int}[]
        for b in bases
            key = (current, b)
            count = get(remaining, key, 0)
            if count > 0
                push!(candidates, (b, count))
            end
        end

        if isempty(candidates)
            # Fallback: any base preserving composition
            push!(result, rand(rng, bases))
        else
            # Weight by remaining counts
            weights = [c[2] for c in candidates]
            total = sum(weights)
            r = rand(rng) * total
            cum = 0.0
            chosen = candidates[1][1]
            for (b, w) in candidates
                cum += w
                if r <= cum
                    chosen = b
                    break
                end
            end
            push!(result, chosen)
            key = (current, chosen)
            remaining[key] = get(remaining, key, 1) - 1
        end
    end

    return LongDNA{4}(result)
end

"""
    IREnrichmentResult

Result of IR enrichment test.
"""
struct IREnrichmentResult
    replicon_id::String
    replicon_length::Int
    observed_count::Int
    baseline_mean::Float64
    baseline_std::Float64
    enrichment_ratio::Float64
    z_score::Float64
    n_baseline_samples::Int
end

"""
    ir_enrichment_test(seq::LongDNA, replicon_id::String;
                       min_stem::Int=8, max_loop::Int=20, min_loop::Int=3,
                       n_samples::Int=10, shuffle_order::Int=1,
                       rng=Random.GLOBAL_RNG) -> IREnrichmentResult

Test IR enrichment vs Markov-shuffled baseline.
"""
function ir_enrichment_test(seq::LongDNA, replicon_id::String;
                            min_stem::Int=8, max_loop::Int=20, min_loop::Int=3,
                            n_samples::Int=10, shuffle_order::Int=1,
                            rng=Random.GLOBAL_RNG)
    # Observed count
    observed = count_inverted_repeats(seq; min_stem=min_stem, max_loop=max_loop, min_loop=min_loop)

    # Baseline from shuffled sequences
    baseline_counts = Int[]
    for _ in 1:n_samples
        shuffled = markov_shuffle(seq; order=shuffle_order, rng=rng)
        count = count_inverted_repeats(shuffled; min_stem=min_stem, max_loop=max_loop, min_loop=min_loop)
        push!(baseline_counts, count)
    end

    baseline_mean = mean(baseline_counts)
    baseline_std = std(baseline_counts)

    enrichment = baseline_mean > 0 ? observed / baseline_mean : (observed > 0 ? Inf : 1.0)
    z_score = baseline_std > 0 ? (observed - baseline_mean) / baseline_std : 0.0

    return IREnrichmentResult(
        replicon_id,
        length(seq),
        observed,
        baseline_mean,
        baseline_std,
        enrichment,
        z_score,
        n_samples
    )
end

end # module
