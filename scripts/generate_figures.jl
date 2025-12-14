#!/usr/bin/env julia
"""
Generate publication figures from analysis results.

Creates:
- Figure 1: Operator map / group structure diagram
- Figure 2: Symmetry statistics (orbit distributions, fixed points)
"""

using ArgParse
using CSV
using DataFrames
using CairoMakie
using Statistics
using StatsBase

function parse_args()
    s = ArgParseSettings(
        description = "Generate publication figures",
        prog = "generate_figures.jl"
    )

    @add_arg_table! s begin
        "--results", "-r"
            help = "Results directory"
            arg_type = String
            default = "results/"
        "--out", "-o"
            help = "Output directory for figures"
            arg_type = String
            default = "results/figures/"
        "--format", "-f"
            help = "Output format (pdf, png, svg)"
            arg_type = String
            default = "pdf"
    end

    return ArgParse.parse_args(s)
end

"""
Generate Figure 1: Operator algebra structure diagram.
"""
function generate_fig1_operator_map(out_dir::String, format::String)
    fig = Figure(size=(800, 600), fontsize=14)

    # Main title
    Label(fig[0, 1:2], "Genomic Operator Algebra", fontsize=20, font=:bold)

    # Left panel: Group structure (invertible operators)
    ax1 = Axis(fig[1, 1],
        title="Dihedral Group D_n (Invertible)",
        aspect=1,
        limits=(-1.5, 1.5, -1.5, 1.5))
    hidedecorations!(ax1)

    # Draw group elements as nodes
    group_ops = ["id", "S", "S²", "...", "R", "RS", "RS²", "..."]
    angles = range(0, 2π, length=length(group_ops)+1)[1:end-1]
    r = 1.0

    for (i, (op, θ)) in enumerate(zip(group_ops, angles))
        x, y = r * cos(θ), r * sin(θ)
        scatter!(ax1, [x], [y], markersize=30, color=:steelblue)
        text!(ax1, x, y, text=op, align=(:center, :center), fontsize=10, color=:white)
    end

    # Center label
    text!(ax1, 0, 0, text="Group\n(closed)", align=(:center, :center), fontsize=12)

    # Right panel: Semigroup with indels
    ax2 = Axis(fig[1, 2],
        title="Semigroup (With Indels)",
        aspect=1,
        limits=(-1.5, 1.5, -1.5, 1.5))
    hidedecorations!(ax2)

    # Draw semigroup elements
    semi_ops = ["S", "R", "K", "M", "D", "I", "V", "..."]
    for (i, (op, θ)) in enumerate(zip(semi_ops, angles))
        x, y = r * cos(θ), r * sin(θ)
        color = op in ["D", "I"] ? :crimson : :steelblue
        scatter!(ax2, [x], [y], markersize=30, color=color)
        text!(ax2, x, y, text=op, align=(:center, :center), fontsize=10, color=:white)
    end

    # Center label
    text!(ax2, 0, 0, text="Semigroup\n(not closed)", align=(:center, :center), fontsize=12)

    # Legend
    Legend(fig[2, 1:2],
        [MarkerElement(marker=:circle, color=:steelblue, markersize=15),
         MarkerElement(marker=:circle, color=:crimson, markersize=15)],
        ["Invertible operators (S, R, K, M, V)", "Non-invertible operators (D, I)"],
        orientation=:horizontal, framevisible=false)

    # Bottom panel: Key identities
    ax3 = Axis(fig[3, 1:2], limits=(0, 1, 0, 1))
    hidedecorations!(ax3)
    hidespines!(ax3)

    identities_text = """
    Key Identities:
    • R∘R = id    (reflection is involution)
    • K∘K = id    (complement is involution)
    • S^n = id    (n rotations = identity)
    • R∘S = S⁻¹∘R (dihedral relation)
    • D∘I ≠ id    (indels non-invertible)
    """

    text!(ax3, 0.5, 0.5, text=identities_text, align=(:center, :center), fontsize=11)

    filepath = joinpath(out_dir, "fig1_operator_map.$format")
    save(filepath, fig)
    println("Saved: $filepath")

    return fig
end

"""
Generate Figure 2: Symmetry statistics from genome analysis.
"""
function generate_fig2_symmetry_stats(results_dir::String, out_dir::String, format::String)
    # Load data
    window_stats_path = joinpath(results_dir, "tables", "window_stats.csv")
    replicon_stats_path = joinpath(results_dir, "tables", "replicon_stats.csv")
    approx_stats_path = joinpath(results_dir, "tables", "approx_symmetry_stats.csv")

    if !isfile(window_stats_path)
        @warn "Window stats not found: $window_stats_path"
        return nothing
    end

    window_df = CSV.read(window_stats_path, DataFrame)
    replicon_df = isfile(replicon_stats_path) ? CSV.read(replicon_stats_path, DataFrame) : nothing
    approx_df = isfile(approx_stats_path) ? CSV.read(approx_stats_path, DataFrame) : nothing

    fig = Figure(size=(1200, 1000), fontsize=12)

    Label(fig[0, 1:3], "Dihedral Symmetry Analysis in Bacterial Genomes", fontsize=18, font=:bold)

    # Panel A: Orbit ratio distribution
    ax1 = Axis(fig[1, 1],
        xlabel="Orbit Ratio (orbit_size / 2n)",
        ylabel="Density",
        title="(A) Orbit Size Distribution")

    if nrow(window_df) > 0
        hist!(ax1, window_df.orbit_ratio, bins=50, normalization=:pdf, color=(:steelblue, 0.7))
        vlines!(ax1, [1.0], color=:red, linestyle=:dash, linewidth=2, label="Max (no symmetry)")

        mean_ratio = mean(window_df.orbit_ratio)
        vlines!(ax1, [mean_ratio], color=:green, linestyle=:dot, linewidth=2, label="Mean")

        axislegend(ax1, position=:lt)
    end

    # Panel B: Orbit ratio by replicon type
    ax2 = Axis(fig[1, 2],
        xlabel="Replicon Type",
        ylabel="Orbit Ratio",
        title="(B) Symmetry by Replicon Type")

    if nrow(window_df) > 0 && "replicon_type" in names(window_df)
        types = unique(window_df.replicon_type)
        type_colors = Dict("chromosome" => :steelblue, "plasmid" => :orange, "unknown" => :gray)

        for (i, t) in enumerate(types)
            subset = filter(row -> row.replicon_type == t, window_df)
            if nrow(subset) > 0
                boxplot!(ax2, fill(i, nrow(subset)), subset.orbit_ratio,
                    color=get(type_colors, t, :gray), width=0.6, label=t)
            end
        end

        ax2.xticks = (1:length(types), types)
    end

    # Panel C: Exact fixed point incidence (expected ~0 under null)
    ax3 = Axis(fig[2, 1],
        xlabel="Exact Invariance Type",
        ylabel="Fraction of Windows",
        title="(C) Exact Fixed Points (Expected ≈ 0)")

    if nrow(window_df) > 0
        pal_frac = mean(window_df.is_palindrome)
        rc_frac = mean(window_df.is_rc_fixed)

        barplot!(ax3, [1, 2], [pal_frac, rc_frac],
            color=[:steelblue, :orange],
            bar_labels=:y,
            label_formatter=x -> "$(round(x*100, digits=1))%")

        ax3.xticks = ([1, 2], ["Palindrome\n(R-fixed)", "RC-fixed"])
        ylims!(ax3, 0, max(0.1, max(pal_frac, rc_frac) * 1.3))
    end

    # Panel D: Replicon length vs symmetry
    ax4 = Axis(fig[2, 2],
        xlabel="Replicon Length (Mbp)",
        ylabel="Mean Orbit Ratio",
        title="(D) Length vs Symmetry",
        xscale=log10)

    if replicon_df !== nothing && nrow(replicon_df) > 0
        type_colors = Dict("chromosome" => :steelblue, "plasmid" => :orange, "unknown" => :gray)

        for t in unique(replicon_df.replicon_type)
            subset = filter(row -> row.replicon_type == t, replicon_df)
            if nrow(subset) > 0
                scatter!(ax4,
                    subset.replicon_length ./ 1e6,
                    subset.mean_orbit_ratio,
                    color=get(type_colors, t, :gray),
                    markersize=6,
                    alpha=0.6,
                    label=t)
            end
        end

        axislegend(ax4, position=:rb)
    end

    # Panel E: Approximate symmetry (d_min/L) across window sizes
    ax5 = Axis(fig[3, 1:2],
        xlabel="Window Size (bp)",
        ylabel="d_min / L (Normalized Distance)",
        title="(E) Approximate Dihedral Self-Similarity: Real vs Shuffled Baseline")

    if approx_df !== nothing && nrow(approx_df) > 0
        window_sizes = sort(unique(approx_df.window_size))

        # Compute means and stds for each window size
        real_means = Float64[]
        real_stds = Float64[]
        shuf_means = Float64[]
        shuf_stds = Float64[]

        for L in window_sizes
            subset = filter(row -> row.window_size == L, approx_df)
            push!(real_means, mean(subset.d_min_normalized_real))
            push!(real_stds, std(subset.d_min_normalized_real))
            push!(shuf_means, mean(subset.d_min_normalized_shuffled))
            push!(shuf_stds, std(subset.d_min_normalized_shuffled))
        end

        # Plot with error bands
        x_pos = 1:length(window_sizes)

        # Real sequences
        band!(ax5, x_pos, real_means .- real_stds, real_means .+ real_stds,
            color=(:steelblue, 0.3))
        scatterlines!(ax5, x_pos, real_means,
            color=:steelblue, linewidth=2, markersize=10, label="Real genomes")

        # Shuffled baseline
        band!(ax5, x_pos, shuf_means .- shuf_stds, shuf_means .+ shuf_stds,
            color=(:orange, 0.3))
        scatterlines!(ax5, x_pos, shuf_means,
            color=:orange, linewidth=2, markersize=10, linestyle=:dash, label="GC-shuffled")

        ax5.xticks = (x_pos, string.(window_sizes))
        axislegend(ax5, position=:rt)
    else
        text!(ax5, 0.5, 0.5, text="Run analyze_approx_symmetry.jl first",
            align=(:center, :center), fontsize=14)
    end

    filepath = joinpath(out_dir, "fig2_symmetry_stats.$format")
    save(filepath, fig)
    println("Saved: $filepath")

    return fig
end

function main()
    args = parse_args()

    results_dir = args["results"]
    out_dir = args["out"]
    format = args["format"]

    mkpath(out_dir)

    println("\n=== Generating Figures ===")

    # Figure 1: Operator map
    println("\nGenerating Figure 1: Operator map...")
    generate_fig1_operator_map(out_dir, format)

    # Figure 2: Symmetry statistics
    println("\nGenerating Figure 2: Symmetry statistics...")
    generate_fig2_symmetry_stats(results_dir, out_dir, format)

    println("\n=== Figures Complete ===")
end

main()
