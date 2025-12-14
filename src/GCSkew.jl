"""
GC Skew Module for Replication Geometry Estimation

Estimates ori/ter positions via cumulative GC skew analysis.
ori ≈ argmin(cumulative skew), ter ≈ argmax(cumulative skew).

Reference: Lobry (1996), Genome Biology 2000.
"""
module GCSkew

using BioSequences
using Statistics

export
    compute_gc_skew, cumulative_gc_skew, estimate_ori_ter,
    OriTerEstimate, skew_confidence_index, split_replichores

"""
    OriTerEstimate

Estimated replication origin and terminus positions.
"""
struct OriTerEstimate
    ori_pos::Int              # 1-indexed position of estimated ori
    ter_pos::Int              # 1-indexed position of estimated ter
    skew_amplitude::Float64   # Max - min of cumulative skew (normalized)
    confidence::Symbol        # :high, :medium, :low
    cumulative_skew::Vector{Float64}
end

"""
    compute_gc_skew(seq::LongDNA; window_size::Int=1000, step::Int=500) -> Vector{Float64}

Compute windowed GC skew: (G-C)/(G+C) for each window.
Returns vector of skew values for each window position.
"""
function compute_gc_skew(seq::LongDNA; window_size::Int=1000, step::Int=500)
    n = length(seq)
    skew = Float64[]

    pos = 1
    while pos + window_size - 1 <= n
        window = seq[pos:pos+window_size-1]

        g_count = 0
        c_count = 0
        for i in 1:window_size
            b = window[i]
            if b == DNA_G
                g_count += 1
            elseif b == DNA_C
                c_count += 1
            end
        end

        total = g_count + c_count
        if total > 0
            push!(skew, (g_count - c_count) / total)
        else
            push!(skew, 0.0)
        end

        pos += step
    end

    return skew
end

"""
    cumulative_gc_skew(skew::Vector{Float64}) -> Vector{Float64}

Compute cumulative sum of GC skew values.
"""
function cumulative_gc_skew(skew::Vector{Float64})
    return cumsum(skew)
end

"""
    skew_confidence_index(cum_skew::Vector{Float64}, seq_length::Int) -> Tuple{Float64, Symbol}

Compute confidence index for ori/ter estimation.
Based on amplitude of cumulative skew relative to genome size.

Returns (amplitude, confidence_level).
"""
function skew_confidence_index(cum_skew::Vector{Float64}, seq_length::Int)
    isempty(cum_skew) && return (0.0, :low)

    amplitude = maximum(cum_skew) - minimum(cum_skew)
    # Normalize by number of windows
    normalized_amp = amplitude / length(cum_skew)

    # Thresholds based on typical bacterial genomes
    if normalized_amp > 0.3
        return (normalized_amp, :high)
    elseif normalized_amp > 0.1
        return (normalized_amp, :medium)
    else
        return (normalized_amp, :low)
    end
end

"""
    estimate_ori_ter(seq::LongDNA; window_size::Int=1000, step::Int=500) -> OriTerEstimate

Estimate replication origin and terminus positions using cumulative GC skew.
ori ≈ position of minimum cumulative skew
ter ≈ position of maximum cumulative skew
"""
function estimate_ori_ter(seq::LongDNA; window_size::Int=1000, step::Int=500)
    n = length(seq)

    # Compute GC skew
    skew = compute_gc_skew(seq; window_size=window_size, step=step)
    isempty(skew) && return OriTerEstimate(1, div(n, 2), 0.0, :low, Float64[])

    # Cumulative skew
    cum_skew = cumulative_gc_skew(skew)

    # Find min/max positions
    min_idx = argmin(cum_skew)
    max_idx = argmax(cum_skew)

    # Convert window index to sequence position
    ori_pos = (min_idx - 1) * step + div(window_size, 2)
    ter_pos = (max_idx - 1) * step + div(window_size, 2)

    # Clamp to valid range
    ori_pos = clamp(ori_pos, 1, n)
    ter_pos = clamp(ter_pos, 1, n)

    # Confidence
    amplitude, confidence = skew_confidence_index(cum_skew, n)

    return OriTerEstimate(ori_pos, ter_pos, amplitude, confidence, cum_skew)
end

"""
    split_replichores(seq::LongDNA, ori_pos::Int) -> Tuple{LongDNA{4}, LongDNA{4}}

Split circular genome at ori into two replichores.
Returns (leading_replichore, lagging_replichore).

Leading: ori -> ter (first half after rotation)
Lagging: ter -> ori (second half after rotation)
"""
function split_replichores(seq::LongDNA, ori_pos::Int)
    n = length(seq)

    # Rotate so ori is at position 1
    if ori_pos > 1
        rotated = seq[ori_pos:end] * seq[1:ori_pos-1]
    else
        rotated = seq
    end

    half = div(n, 2)
    leading = rotated[1:half]
    lagging = rotated[half+1:end]

    return (leading, lagging)
end

"""
    compute_replichore_gc_content(seq::LongDNA, ori_pos::Int) -> NamedTuple

Compute GC content for each replichore.
"""
function compute_replichore_gc_content(seq::LongDNA, ori_pos::Int)
    leading, lagging = split_replichores(seq, ori_pos)

    function gc_content(s::LongDNA)
        gc = 0
        total = 0
        for i in 1:length(s)
            b = s[i]
            if b == DNA_G || b == DNA_C
                gc += 1
                total += 1
            elseif b == DNA_A || b == DNA_T
                total += 1
            end
        end
        return total > 0 ? gc / total : 0.0
    end

    return (
        leading_gc = gc_content(leading),
        lagging_gc = gc_content(lagging),
        gc_diff = gc_content(leading) - gc_content(lagging)
    )
end

end # module
