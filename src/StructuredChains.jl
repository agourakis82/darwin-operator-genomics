"""
Structured Operator Chains Module

Generates biologically-constrained operator chains mimicking known bacterial
evolutionary patterns: symmetric inversions around ori/ter (X-alignment mechanism).

Reference: Genome Biology 2000 (X-alignments).
"""
module StructuredChains

using Random

export
    StructuredChainGenerator, generate_ori_symmetric_chain,
    generate_ter_symmetric_chain, generate_mixed_structured_chain,
    chain_to_symbols, OperatorEvent

"""
    OperatorEvent

Represents a single operator in a chain.
"""
struct OperatorEvent
    op_type::Symbol           # :shift, :reverse, :inversion, :mutation
    param1::Int               # Position or shift amount
    param2::Int               # End position for inversion, or 0
    symmetric::Bool           # Whether this is symmetric around ori/ter
end

"""
    StructuredChainGenerator

Generator for structured operator chains.
"""
struct StructuredChainGenerator
    genome_length::Int
    ori_pos::Int
    ter_pos::Int
    rng::AbstractRNG
end

function StructuredChainGenerator(genome_length::Int, ori_pos::Int, ter_pos::Int; seed::Int=42)
    return StructuredChainGenerator(genome_length, ori_pos, ter_pos, MersenneTwister(seed))
end

"""
    generate_ori_symmetric_chain(gen::StructuredChainGenerator, n_events::Int;
                                  inversion_prob::Float64=0.7,
                                  mutation_rate::Float64=0.1) -> Vector{OperatorEvent}

Generate a chain with symmetric inversions around ori.
This mimics the X-alignment pattern observed in bacterial evolution.
"""
function generate_ori_symmetric_chain(gen::StructuredChainGenerator, n_events::Int;
                                       inversion_prob::Float64=0.7,
                                       mutation_rate::Float64=0.1)
    events = OperatorEvent[]
    n = gen.genome_length
    ori = gen.ori_pos

    for _ in 1:n_events
        if rand(gen.rng) < inversion_prob
            # Symmetric inversion around ori
            # Choose distance from ori
            max_dist = min(ori - 1, n - ori)
            max_dist < 100 && continue

            dist = rand(gen.rng, 100:max_dist)
            inv_len = rand(gen.rng, 50:min(500, dist))

            # Positions symmetric around ori
            pos1 = ori - dist
            pos2 = ori + dist

            # Invert both segments (X-alignment style)
            push!(events, OperatorEvent(:inversion, pos1, pos1 + inv_len, true))
            push!(events, OperatorEvent(:inversion, pos2 - inv_len, pos2, true))
        else
            # Random small mutation
            if rand(gen.rng) < mutation_rate
                pos = rand(gen.rng, 1:n)
                push!(events, OperatorEvent(:mutation, pos, 0, false))
            else
                # Small shift
                shift = rand(gen.rng, 1:10)
                push!(events, OperatorEvent(:shift, shift, 0, false))
            end
        end
    end

    return events
end

"""
    generate_ter_symmetric_chain(gen::StructuredChainGenerator, n_events::Int;
                                  inversion_prob::Float64=0.7,
                                  mutation_rate::Float64=0.1) -> Vector{OperatorEvent}

Generate a chain with symmetric inversions around ter.
"""
function generate_ter_symmetric_chain(gen::StructuredChainGenerator, n_events::Int;
                                       inversion_prob::Float64=0.7,
                                       mutation_rate::Float64=0.1)
    events = OperatorEvent[]
    n = gen.genome_length
    ter = gen.ter_pos

    for _ in 1:n_events
        if rand(gen.rng) < inversion_prob
            # Symmetric inversion around ter
            max_dist = min(ter - 1, n - ter)
            max_dist < 100 && continue

            dist = rand(gen.rng, 100:max_dist)
            inv_len = rand(gen.rng, 50:min(500, dist))

            pos1 = ter - dist
            pos2 = ter + dist

            push!(events, OperatorEvent(:inversion, pos1, pos1 + inv_len, true))
            push!(events, OperatorEvent(:inversion, pos2 - inv_len, pos2, true))
        else
            if rand(gen.rng) < mutation_rate
                pos = rand(gen.rng, 1:n)
                push!(events, OperatorEvent(:mutation, pos, 0, false))
            else
                shift = rand(gen.rng, 1:10)
                push!(events, OperatorEvent(:shift, shift, 0, false))
            end
        end
    end

    return events
end

"""
    generate_mixed_structured_chain(gen::StructuredChainGenerator, n_events::Int;
                                     ori_weight::Float64=0.5) -> Vector{OperatorEvent}

Generate chain with symmetric inversions around both ori and ter.
"""
function generate_mixed_structured_chain(gen::StructuredChainGenerator, n_events::Int;
                                          ori_weight::Float64=0.5)
    events = OperatorEvent[]

    for _ in 1:n_events
        if rand(gen.rng) < ori_weight
            chain = generate_ori_symmetric_chain(gen, 1)
        else
            chain = generate_ter_symmetric_chain(gen, 1)
        end
        append!(events, chain)
    end

    return events
end

"""
    chain_to_symbols(events::Vector{OperatorEvent}) -> Vector{Symbol}

Convert event chain to symbol sequence for comparison with existing analysis.
"""
function chain_to_symbols(events::Vector{OperatorEvent})
    symbols = Symbol[]
    for e in events
        push!(symbols, e.op_type)
    end
    return symbols
end

"""
    count_symmetric_events(events::Vector{OperatorEvent}) -> Int

Count number of symmetric (ori/ter-balanced) events in chain.
"""
function count_symmetric_events(events::Vector{OperatorEvent})
    return count(e -> e.symmetric, events)
end

"""
    chain_complexity(events::Vector{OperatorEvent}) -> NamedTuple

Compute complexity metrics for a chain.
"""
function chain_complexity(events::Vector{OperatorEvent})
    n = length(events)
    n_symmetric = count_symmetric_events(events)
    n_inversions = count(e -> e.op_type == :inversion, events)
    n_mutations = count(e -> e.op_type == :mutation, events)
    n_shifts = count(e -> e.op_type == :shift, events)

    return (
        total_events = n,
        symmetric_events = n_symmetric,
        symmetric_ratio = n > 0 ? n_symmetric / n : 0.0,
        inversions = n_inversions,
        mutations = n_mutations,
        shifts = n_shifts
    )
end

end # module
