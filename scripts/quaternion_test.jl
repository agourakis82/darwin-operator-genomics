#!/usr/bin/env julia
"""
Quaternionic compression hypothesis test.

Tests whether quaternion-based representations of operator chains
outperform Markov baselines for next-operator prediction.
"""

using ArgParse
using JSON3
using FASTX
using BioSequences
using CodecZlib
using Random
using Statistics
using StatsBase
using DataFrames
using CSV
using CairoMakie
using Rotations
using ProgressMeter
using Dates

using DarwinOperatorGenomics

function parse_args()
    s = ArgParseSettings(
        description = "Test quaternionic compression hypothesis",
        prog = "quaternion_test.jl"
    )

    @add_arg_table! s begin
        "--cache", "-c"
            help = "Cache directory with downloaded genomes"
            arg_type = String
            default = "data/cache/"
        "--out", "-o"
            help = "Output directory for results"
            arg_type = String
            default = "results/"
        "--n-chains"
            help = "Number of operator chains to generate"
            arg_type = Int
            default = 1000
        "--chain-length"
            help = "Length of each operator chain"
            arg_type = Int
            default = 50
        "--seed", "-s"
            help = "Random seed"
            arg_type = Int
            default = 42
    end

    return ArgParse.parse_args(s)
end

# =============================================================================
# Baseline Models
# =============================================================================

"""
Markov model for operator sequences.
Predicts next operator based on previous k operators.
"""
struct MarkovModel
    order::Int
    transition_counts::Dict{Vector{Int}, Vector{Int}}
    total_counts::Dict{Vector{Int}, Int}
end

function MarkovModel(order::Int)
    MarkovModel(order, Dict{Vector{Int}, Vector{Int}}(), Dict{Vector{Int}, Int}())
end

function train!(model::MarkovModel, sequences::Vector{Vector{Int}}, n_tokens::Int)
    for seq in sequences
        for i in (model.order + 1):length(seq)
            context = seq[i-model.order:i-1]
            target = seq[i]

            if !haskey(model.transition_counts, context)
                model.transition_counts[context] = zeros(Int, n_tokens)
                model.total_counts[context] = 0
            end

            model.transition_counts[context][target] += 1
            model.total_counts[context] += 1
        end
    end
end

function predict_probs(model::MarkovModel, context::Vector{Int}, n_tokens::Int)
    if haskey(model.transition_counts, context) && model.total_counts[context] > 0
        counts = model.transition_counts[context]
        total = model.total_counts[context]
        # Add smoothing
        return (counts .+ 1) ./ (total + n_tokens)
    else
        # Uniform distribution
        return fill(1.0 / n_tokens, n_tokens)
    end
end

function evaluate_markov(model::MarkovModel, sequences::Vector{Vector{Int}}, n_tokens::Int)
    correct = 0
    total = 0
    log_likelihood = 0.0

    for seq in sequences
        for i in (model.order + 1):length(seq)
            context = seq[i-model.order:i-1]
            target = seq[i]

            probs = predict_probs(model, context, n_tokens)
            predicted = argmax(probs)

            if predicted == target
                correct += 1
            end

            log_likelihood += log(probs[target] + 1e-10)
            total += 1
        end
    end

    accuracy = correct / total
    perplexity = exp(-log_likelihood / total)

    return (accuracy=accuracy, perplexity=perplexity)
end

# =============================================================================
# Quaternion Model
# =============================================================================

"""
Quaternion-based operator model.
Maintains a latent quaternion state that gets updated by operator embeddings.
"""
struct QuaternionModel
    operator_quats::Vector{QuatRotation{Float64}}  # Embedding for each operator
    readout_weights::Matrix{Float64}  # Maps quaternion to prediction logits
end

function QuaternionModel(n_tokens::Int; rng=Random.GLOBAL_RNG)
    # Initialize operator quaternions as small random rotations
    operator_quats = [QuatRotation(randn(rng, 4)...) for _ in 1:n_tokens]

    # Readout: 4D quaternion -> n_tokens logits
    readout_weights = randn(rng, n_tokens, 4) * 0.1

    return QuaternionModel(operator_quats, readout_weights)
end

function forward_quaternion(model::QuaternionModel, sequence::Vector{Int})
    # Start with identity quaternion
    state = QuatRotation(1.0, 0.0, 0.0, 0.0)

    states = QuatRotation{Float64}[]
    push!(states, state)

    for op_idx in sequence[1:end-1]
        # Update state by quaternion multiplication
        op_quat = model.operator_quats[op_idx]
        state = state * op_quat
        push!(states, state)
    end

    return states
end

function predict_quaternion(model::QuaternionModel, state::QuatRotation)
    # Extract quaternion components
    q = [state.q.s, state.q.v1, state.q.v2, state.q.v3]

    # Linear readout
    logits = model.readout_weights * q

    # Softmax
    exp_logits = exp.(logits .- maximum(logits))
    probs = exp_logits ./ sum(exp_logits)

    return probs
end

function evaluate_quaternion(model::QuaternionModel, sequences::Vector{Vector{Int}}, n_tokens::Int)
    correct = 0
    total = 0
    log_likelihood = 0.0

    for seq in sequences
        states = forward_quaternion(model, seq)

        for i in 2:length(seq)
            state = states[i]
            target = seq[i]

            probs = predict_quaternion(model, state)
            predicted = argmax(probs)

            if predicted == target
                correct += 1
            end

            log_likelihood += log(probs[target] + 1e-10)
            total += 1
        end
    end

    accuracy = correct / total
    perplexity = exp(-log_likelihood / total)

    return (accuracy=accuracy, perplexity=perplexity)
end

# Simple gradient-free optimization for quaternion model
function train_quaternion!(model::QuaternionModel, train_seqs::Vector{Vector{Int}},
                          n_tokens::Int; n_iters=100, lr=0.01, rng=Random.GLOBAL_RNG)
    best_accuracy = 0.0

    for iter in 1:n_iters
        # Random perturbation
        for i in 1:n_tokens
            delta = randn(rng, 4) * lr / sqrt(iter)
            old_q = model.operator_quats[i]
            new_q_params = [old_q.q.s, old_q.q.v1, old_q.q.v2, old_q.q.v3] .+ delta
            model.operator_quats[i] = QuatRotation(new_q_params...)
        end

        # Perturb readout
        model.readout_weights .+= randn(rng, size(model.readout_weights)...) * lr / sqrt(iter)

        result = evaluate_quaternion(model, train_seqs, n_tokens)
        if result.accuracy > best_accuracy
            best_accuracy = result.accuracy
        end
    end

    return best_accuracy
end

# =============================================================================
# Main Experiment
# =============================================================================

function load_manifest(cache_dir::String)
    manifest_path = joinpath(cache_dir, "manifest.jsonl")
    entries = Dict[]

    if isfile(manifest_path)
        for line in eachline(manifest_path)
            isempty(strip(line)) && continue
            push!(entries, JSON3.read(line, Dict))
        end
    end

    return entries
end

function generate_operator_sequences(n_chains::Int, chain_length::Int, seq_len::Int;
                                     rng=Random.GLOBAL_RNG, include_indels::Bool=false)
    sequences = Vector{Int}[]

    for _ in 1:n_chains
        chain = generate_random_chain(seq_len, chain_length; rng=rng, include_indels=include_indels)
        tokens = [operator_to_token(op) for op in chain]
        push!(sequences, tokens)
    end

    return sequences
end

function main()
    args = parse_args()

    cache_dir = args["cache"]
    out_dir = args["out"]
    n_chains = args["n-chains"]
    chain_length = args["chain-length"]
    seed = args["seed"]

    rng = MersenneTwister(seed)

    mkpath(out_dir)
    mkpath(joinpath(out_dir, "figures"))
    mkpath(joinpath(out_dir, "text"))

    println("\n=== Quaternionic Compression Hypothesis Test ===")

    # Generate operator sequences
    println("\nGenerating operator chains...")
    println("  N chains: $n_chains")
    println("  Chain length: $chain_length")

    # Without indels (group case)
    group_seqs = generate_operator_sequences(n_chains, chain_length, 1000;
                                              rng=rng, include_indels=false)

    # With indels (semigroup case)
    semi_seqs = generate_operator_sequences(n_chains, chain_length, 1000;
                                             rng=rng, include_indels=true)

    # Split into train/test
    n_train = div(n_chains * 4, 5)

    group_train = group_seqs[1:n_train]
    group_test = group_seqs[n_train+1:end]

    semi_train = semi_seqs[1:n_train]
    semi_test = semi_seqs[n_train+1:end]

    n_tokens_group = 5  # S, R, K, RC, M (without indels)
    n_tokens_semi = NUM_OPERATOR_TYPES  # All operators

    # Results storage
    results = DataFrame(
        condition = String[],
        model = String[],
        accuracy = Float64[],
        perplexity = Float64[]
    )

    # ==========================================================================
    # Experiment 1: Group case (invertible operators only)
    # ==========================================================================
    println("\n--- Group Case (Invertible Operators) ---")

    # Markov order 1
    println("  Training Markov(1)...")
    markov1_group = MarkovModel(1)
    train!(markov1_group, group_train, n_tokens_group)
    m1g_result = evaluate_markov(markov1_group, group_test, n_tokens_group)
    push!(results, ("group", "Markov(1)", m1g_result.accuracy, m1g_result.perplexity))
    println("    Accuracy: $(round(m1g_result.accuracy * 100, digits=1))%")

    # Markov order 2
    println("  Training Markov(2)...")
    markov2_group = MarkovModel(2)
    train!(markov2_group, group_train, n_tokens_group)
    m2g_result = evaluate_markov(markov2_group, group_test, n_tokens_group)
    push!(results, ("group", "Markov(2)", m2g_result.accuracy, m2g_result.perplexity))
    println("    Accuracy: $(round(m2g_result.accuracy * 100, digits=1))%")

    # Quaternion model
    println("  Training Quaternion...")
    quat_group = QuaternionModel(n_tokens_group; rng=rng)
    train_quaternion!(quat_group, group_train, n_tokens_group; n_iters=50, rng=rng)
    qg_result = evaluate_quaternion(quat_group, group_test, n_tokens_group)
    push!(results, ("group", "Quaternion", qg_result.accuracy, qg_result.perplexity))
    println("    Accuracy: $(round(qg_result.accuracy * 100, digits=1))%")

    # ==========================================================================
    # Experiment 2: Semigroup case (with indels)
    # ==========================================================================
    println("\n--- Semigroup Case (With Indels) ---")

    # Markov order 1
    println("  Training Markov(1)...")
    markov1_semi = MarkovModel(1)
    train!(markov1_semi, semi_train, n_tokens_semi)
    m1s_result = evaluate_markov(markov1_semi, semi_test, n_tokens_semi)
    push!(results, ("semigroup", "Markov(1)", m1s_result.accuracy, m1s_result.perplexity))
    println("    Accuracy: $(round(m1s_result.accuracy * 100, digits=1))%")

    # Markov order 2
    println("  Training Markov(2)...")
    markov2_semi = MarkovModel(2)
    train!(markov2_semi, semi_train, n_tokens_semi)
    m2s_result = evaluate_markov(markov2_semi, semi_test, n_tokens_semi)
    push!(results, ("semigroup", "Markov(2)", m2s_result.accuracy, m2s_result.perplexity))
    println("    Accuracy: $(round(m2s_result.accuracy * 100, digits=1))%")

    # Quaternion model
    println("  Training Quaternion...")
    quat_semi = QuaternionModel(n_tokens_semi; rng=rng)
    train_quaternion!(quat_semi, semi_train, n_tokens_semi; n_iters=50, rng=rng)
    qs_result = evaluate_quaternion(quat_semi, semi_test, n_tokens_semi)
    push!(results, ("semigroup", "Quaternion", qs_result.accuracy, qs_result.perplexity))
    println("    Accuracy: $(round(qs_result.accuracy * 100, digits=1))%")

    # ==========================================================================
    # Save results
    # ==========================================================================
    println("\n=== Saving Results ===")

    results_path = joinpath(out_dir, "tables", "quaternion_results.csv")
    CSV.write(results_path, results)
    println("Saved: $results_path")

    # Generate Figure 3
    println("\nGenerating Figure 3...")

    fig = Figure(size=(900, 400), fontsize=12)

    Label(fig[0, 1:2], "Quaternionic Compression Hypothesis Test", fontsize=16, font=:bold)

    # Panel A: Group case
    ax1 = Axis(fig[1, 1],
        xlabel="Model",
        ylabel="Accuracy (%)",
        title="(A) Group Case (Invertible Ops)",
        xticks=(1:3, ["Markov(1)", "Markov(2)", "Quaternion"]))

    group_results = filter(row -> row.condition == "group", results)
    barplot!(ax1, 1:3, group_results.accuracy .* 100,
        color=[:steelblue, :steelblue, :orange])
    ylims!(ax1, 0, 100)

    # Panel B: Semigroup case
    ax2 = Axis(fig[1, 2],
        xlabel="Model",
        ylabel="Accuracy (%)",
        title="(B) Semigroup Case (With Indels)",
        xticks=(1:3, ["Markov(1)", "Markov(2)", "Quaternion"]))

    semi_results = filter(row -> row.condition == "semigroup", results)
    barplot!(ax2, 1:3, semi_results.accuracy .* 100,
        color=[:steelblue, :steelblue, :orange])
    ylims!(ax2, 0, 100)

    Legend(fig[2, 1:2],
        [PolyElement(color=:steelblue), PolyElement(color=:orange)],
        ["Markov Baseline", "Quaternion Model"],
        orientation=:horizontal, framevisible=false)

    fig_path = joinpath(out_dir, "figures", "fig3_quaternion_vs_baselines.pdf")
    save(fig_path, fig)
    println("Saved: $fig_path")

    # Generate text summary
    text_summary = """
# Results: Quaternionic Compression Hypothesis

## Experimental Setup

- **N chains**: $n_chains ($(n_train) train, $(n_chains - n_train) test)
- **Chain length**: $chain_length operators
- **Seed**: $seed

## Results

### Group Case (Invertible Operators Only: S, R, K, RC, M)

| Model | Accuracy | Perplexity |
|-------|----------|------------|
| Markov(1) | $(round(m1g_result.accuracy * 100, digits=1))% | $(round(m1g_result.perplexity, digits=2)) |
| Markov(2) | $(round(m2g_result.accuracy * 100, digits=1))% | $(round(m2g_result.perplexity, digits=2)) |
| Quaternion | $(round(qg_result.accuracy * 100, digits=1))% | $(round(qg_result.perplexity, digits=2)) |

### Semigroup Case (With Indels: D, I, V)

| Model | Accuracy | Perplexity |
|-------|----------|------------|
| Markov(1) | $(round(m1s_result.accuracy * 100, digits=1))% | $(round(m1s_result.perplexity, digits=2)) |
| Markov(2) | $(round(m2s_result.accuracy * 100, digits=1))% | $(round(m2s_result.perplexity, digits=2)) |
| Quaternion | $(round(qs_result.accuracy * 100, digits=1))% | $(round(qs_result.perplexity, digits=2)) |

## Interpretation

"""

    # Add interpretation based on results
    group_quat_better = qg_result.accuracy > max(m1g_result.accuracy, m2g_result.accuracy)
    semi_quat_better = qs_result.accuracy > max(m1s_result.accuracy, m2s_result.accuracy)

    if group_quat_better && semi_quat_better
        text_summary *= """
The quaternion model outperforms Markov baselines in both cases, supporting the
hypothesis that quaternionic representations capture structure in operator sequences.
"""
    elseif group_quat_better
        text_summary *= """
The quaternion model outperforms Markov baselines only in the group case (invertible operators).
This suggests quaternionic structure may be more relevant for algebraically closed operations.
"""
    elseif semi_quat_better
        text_summary *= """
The quaternion model outperforms Markov baselines only in the semigroup case.
This unexpected result warrants further investigation.
"""
    else
        text_summary *= """
**Negative Result**: The quaternion model does not consistently outperform Markov baselines.

This suggests that simple quaternionic representations may not capture meaningful structure
in random operator chains. Possible explanations:

1. The quaternion model lacks sufficient capacity or training
2. Random operator chains lack the sequential structure that quaternions might capture
3. The quaternionic hypothesis may not apply to this domain

This negative result is preserved for scientific completeness. Future work could explore:
- More sophisticated quaternion architectures
- Real evolutionary operator chains (not random)
- Alternative algebraic representations (Clifford algebras, etc.)
"""
    end

    text_summary *= """

---
*Generated: $(now())*
"""

    text_path = joinpath(out_dir, "text", "results_r3.md")
    write(text_path, text_summary)
    println("Saved: $text_path")

    println("\n=== Quaternion Test Complete ===")
end

main()
