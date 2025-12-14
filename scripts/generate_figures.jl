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

"""
Generate Figure 4: K-mer Inversion Symmetry Analysis.
"""
function generate_fig4_kmer_symmetry(results_dir::String, out_dir::String, format::String)
    kmer_path = joinpath(results_dir, "tables", "kmer_inversion_symmetry.csv")
    by_k_path = joinpath(results_dir, "tables", "kmer_inversion_symmetry_by_k.csv")

    if !isfile(kmer_path)
        @warn "K-mer symmetry data not found: $kmer_path"
        return nothing
    end

    kmer_df = CSV.read(kmer_path, DataFrame)
    by_k_df = isfile(by_k_path) ? CSV.read(by_k_path, DataFrame) : nothing

    fig = Figure(size=(1200, 800), fontsize=12)

    Label(fig[0, 1:3], "K-mer Inversion Symmetry (Second Chargaff Parity)", fontsize=16, font=:bold)

    # Panel A: X_k profile (mean across genomes)
    ax1 = Axis(fig[1, 1],
        xlabel="k",
        ylabel="X_k (Inversion Score)",
        title="(A) Mean X_k Profile")

    if by_k_df !== nothing && nrow(by_k_df) > 0
        k_values = sort(unique(by_k_df.k))
        means = [mean(by_k_df[by_k_df.k .== k, :mean_X_k]) for k in k_values]
        stds = [std(by_k_df[by_k_df.k .== k, :mean_X_k]) for k in k_values]

        band!(ax1, k_values, means .- stds, means .+ stds, color=(:steelblue, 0.3))
        scatterlines!(ax1, k_values, means, color=:steelblue, linewidth=2, markersize=8)

        # Tolerance lines
        hlines!(ax1, [0.1], color=:red, linestyle=:dash, label="τ = 0.1")
        hlines!(ax1, [0.05], color=:orange, linestyle=:dot, label="τ = 0.05")
        axislegend(ax1, position=:lt)
    end

    # Panel B: Distribution of symmetry limits
    ax2 = Axis(fig[1, 2],
        xlabel="K_L (Symmetry Limit, τ=0.1)",
        ylabel="Count",
        title="(B) Distribution of K_L")

    if nrow(kmer_df) > 0 && "K_L_tau01" in names(kmer_df)
        hist!(ax2, kmer_df.K_L_tau01, bins=20, color=(:steelblue, 0.7))
    end

    # Panel C: K_L vs genome length
    ax3 = Axis(fig[1, 3],
        xlabel="Replicon Length (Mbp)",
        ylabel="K_L (τ=0.1)",
        title="(C) Symmetry Limit vs Length",
        xscale=log10)

    if nrow(kmer_df) > 0 && "replicon_length" in names(kmer_df) && "K_L_tau01" in names(kmer_df)
        scatter!(ax3, kmer_df.replicon_length ./ 1e6, kmer_df.K_L_tau01,
            color=(:steelblue, 0.5), markersize=6)
    end

    # Panel D: X_1 distribution (base-level Chargaff)
    ax4 = Axis(fig[2, 1],
        xlabel="X_1 (Base-level)",
        ylabel="Count",
        title="(D) Base-Level Symmetry")

    if nrow(kmer_df) > 0 && "X_1" in names(kmer_df)
        hist!(ax4, filter(isfinite, kmer_df.X_1), bins=30, color=(:green, 0.7))
        vlines!(ax4, [0.0], color=:red, linestyle=:dash, label="Perfect symmetry")
        axislegend(ax4, position=:rt)
    end

    # Panel E: X_5 vs X_10
    ax5 = Axis(fig[2, 2],
        xlabel="X_5 (5-mer)",
        ylabel="X_10 (10-mer)",
        title="(E) Multi-Scale Symmetry")

    if nrow(kmer_df) > 0 && "X_5" in names(kmer_df) && "X_10" in names(kmer_df)
        valid = filter(row -> isfinite(row.X_5) && isfinite(row.X_10), kmer_df)
        scatter!(ax5, valid.X_5, valid.X_10, color=(:steelblue, 0.5), markersize=6)
        lines!(ax5, [0, 1], [0, 1], color=:gray, linestyle=:dash, label="y=x")
        axislegend(ax5, position=:rb)
    end

    # Panel F: Summary text
    ax6 = Axis(fig[2, 3], limits=(0, 1, 0, 1))
    hidedecorations!(ax6)
    hidespines!(ax6)

    if nrow(kmer_df) > 0
        median_K_L = "K_L_tau01" in names(kmer_df) ? round(median(kmer_df.K_L_tau01), digits=1) : "N/A"
        summary = """
        Summary Statistics:

        N replicons: $(nrow(kmer_df))
        Median K_L(τ=0.1): $median_K_L

        K_L indicates the largest k where
        k-mer counts satisfy approximate
        Chargaff parity (X_k ≤ τ).
        """
        text!(ax6, 0.5, 0.5, text=summary, align=(:center, :center), fontsize=11)
    end

    filepath = joinpath(out_dir, "fig4_kmer_symmetry.$format")
    save(filepath, fig)
    println("Saved: $filepath")

    return fig
end

"""
Generate Figure 5: GC Skew and Replichore Analysis.
"""
function generate_fig5_gc_skew(results_dir::String, out_dir::String, format::String)
    ori_ter_path = joinpath(results_dir, "tables", "gc_skew_ori_ter.csv")
    replichore_path = joinpath(results_dir, "tables", "replichore_symmetry.csv")

    if !isfile(ori_ter_path)
        @warn "GC skew data not found: $ori_ter_path"
        return nothing
    end

    ori_ter_df = CSV.read(ori_ter_path, DataFrame)
    replichore_df = isfile(replichore_path) ? CSV.read(replichore_path, DataFrame) : nothing

    fig = Figure(size=(1200, 800), fontsize=12)

    Label(fig[0, 1:3], "GC Skew and Replichore Analysis", fontsize=16, font=:bold)

    # Panel A: Confidence distribution
    ax1 = Axis(fig[1, 1],
        xlabel="Confidence Level",
        ylabel="Count",
        title="(A) Ori/Ter Estimation Confidence")

    if nrow(ori_ter_df) > 0 && "confidence" in names(ori_ter_df)
        conf_counts = combine(groupby(ori_ter_df, :confidence), nrow => :count)
        conf_order = ["high", "medium", "low"]
        colors = [:green, :orange, :red]

        positions = Int[]
        counts = Int[]
        bar_colors = Symbol[]

        for (i, c) in enumerate(conf_order)
            row = filter(r -> r.confidence == c, conf_counts)
            if nrow(row) > 0
                push!(positions, i)
                push!(counts, row[1, :count])
                push!(bar_colors, colors[i])
            end
        end

        barplot!(ax1, positions, counts, color=bar_colors)
        ax1.xticks = (1:3, conf_order)
    end

    # Panel B: Skew amplitude distribution
    ax2 = Axis(fig[1, 2],
        xlabel="Skew Amplitude",
        ylabel="Count",
        title="(B) GC Skew Amplitude Distribution")

    if nrow(ori_ter_df) > 0 && "skew_amplitude" in names(ori_ter_df)
        hist!(ax2, filter(isfinite, ori_ter_df.skew_amplitude), bins=30, color=(:steelblue, 0.7))
    end

    # Panel C: Leading vs Lagging GC content
    ax3 = Axis(fig[1, 3],
        xlabel="Leading Strand GC",
        ylabel="Lagging Strand GC",
        title="(C) Replichore GC Content")

    if nrow(ori_ter_df) > 0 && "leading_gc" in names(ori_ter_df)
        scatter!(ax3, ori_ter_df.leading_gc, ori_ter_df.lagging_gc,
            color=(:steelblue, 0.5), markersize=6)
        lines!(ax3, [0.2, 0.8], [0.2, 0.8], color=:gray, linestyle=:dash)
    end

    # Panel D: Replichore X_5 comparison (if available)
    ax4 = Axis(fig[2, 1],
        xlabel="Replichore",
        ylabel="X_5 (Inversion Score)",
        title="(D) Leading vs Lagging Symmetry")

    if replichore_df !== nothing && nrow(replichore_df) > 0 && "X_5" in names(replichore_df)
        leading = filter(r -> r.replichore == "leading", replichore_df)
        lagging = filter(r -> r.replichore == "lagging", replichore_df)

        if nrow(leading) > 0 && nrow(lagging) > 0
            violin!(ax4, fill(1, nrow(leading)), filter(isfinite, leading.X_5), color=(:steelblue, 0.7))
            violin!(ax4, fill(2, nrow(lagging)), filter(isfinite, lagging.X_5), color=(:orange, 0.7))
            ax4.xticks = ([1, 2], ["Leading", "Lagging"])
        end
    end

    # Panel E: GC diff distribution
    ax5 = Axis(fig[2, 2],
        xlabel="GC Diff (Leading - Lagging)",
        ylabel="Count",
        title="(E) GC Content Asymmetry")

    if nrow(ori_ter_df) > 0 && "gc_diff" in names(ori_ter_df)
        hist!(ax5, filter(isfinite, ori_ter_df.gc_diff), bins=30, color=(:purple, 0.7))
        vlines!(ax5, [0], color=:red, linestyle=:dash)
    end

    # Panel F: Summary
    ax6 = Axis(fig[2, 3], limits=(0, 1, 0, 1))
    hidedecorations!(ax6)
    hidespines!(ax6)

    if nrow(ori_ter_df) > 0
        n_high = sum(ori_ter_df.confidence .== "high")
        n_med = sum(ori_ter_df.confidence .== "medium")
        summary = """
        Summary:

        Total replicons: $(nrow(ori_ter_df))
        High confidence: $n_high
        Medium confidence: $n_med

        GC skew estimates replication
        geometry (ori ≈ min, ter ≈ max
        of cumulative skew).
        """
        text!(ax6, 0.5, 0.5, text=summary, align=(:center, :center), fontsize=11)
    end

    filepath = joinpath(out_dir, "fig5_gc_skew.$format")
    save(filepath, fig)
    println("Saved: $filepath")

    return fig
end

"""
Generate Figure 6: Inverted Repeats Enrichment.
"""
function generate_fig6_inverted_repeats(results_dir::String, out_dir::String, format::String)
    ir_path = joinpath(results_dir, "tables", "ir_enrichment_summary.csv")

    if !isfile(ir_path)
        @warn "IR enrichment data not found: $ir_path"
        return nothing
    end

    ir_df = CSV.read(ir_path, DataFrame)

    fig = Figure(size=(1200, 600), fontsize=12)

    Label(fig[0, 1:3], "Inverted Repeats Enrichment Analysis", fontsize=16, font=:bold)

    # Panel A: Z-score distribution
    ax1 = Axis(fig[1, 1],
        xlabel="Z-score",
        ylabel="Count",
        title="(A) IR Enrichment Z-scores")

    if nrow(ir_df) > 0 && "z_score" in names(ir_df)
        valid_z = filter(isfinite, ir_df.z_score)
        if !isempty(valid_z)
            hist!(ax1, valid_z, bins=30, color=(:steelblue, 0.7))
            vlines!(ax1, [2.0], color=:red, linestyle=:dash, label="z = 2")
            vlines!(ax1, [0.0], color=:gray, linestyle=:dot, label="z = 0")
            axislegend(ax1, position=:rt)
        end
    end

    # Panel B: Enrichment ratio distribution
    ax2 = Axis(fig[1, 2],
        xlabel="Enrichment Ratio (Observed / Baseline)",
        ylabel="Count",
        title="(B) IR Enrichment Ratios")

    if nrow(ir_df) > 0 && "enrichment_ratio" in names(ir_df)
        valid_e = filter(x -> isfinite(x) && x < 10, ir_df.enrichment_ratio)
        if !isempty(valid_e)
            hist!(ax2, valid_e, bins=30, color=(:orange, 0.7))
            vlines!(ax2, [1.0], color=:red, linestyle=:dash, label="No enrichment")
            axislegend(ax2, position=:rt)
        end
    end

    # Panel C: Observed vs Expected
    ax3 = Axis(fig[1, 3],
        xlabel="Baseline Mean",
        ylabel="Observed Count",
        title="(C) Observed vs Expected IR Count")

    if nrow(ir_df) > 0 && "observed_ir_count" in names(ir_df) && "baseline_mean" in names(ir_df)
        scatter!(ax3, ir_df.baseline_mean, ir_df.observed_ir_count,
            color=(:steelblue, 0.5), markersize=5)
        # Add y=x line
        max_val = max(maximum(ir_df.baseline_mean), maximum(ir_df.observed_ir_count))
        lines!(ax3, [0, max_val], [0, max_val], color=:gray, linestyle=:dash, label="y=x")
        axislegend(ax3, position=:rb)
    end

    filepath = joinpath(out_dir, "fig6_inverted_repeats.$format")
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

    # Figure 4: K-mer symmetry (Figure 3 is generated by quaternion_test.jl)
    println("\nGenerating Figure 4: K-mer inversion symmetry...")
    generate_fig4_kmer_symmetry(results_dir, out_dir, format)

    # Figure 5: GC skew
    println("\nGenerating Figure 5: GC skew and replichore analysis...")
    generate_fig5_gc_skew(results_dir, out_dir, format)

    # Figure 6: Inverted repeats
    println("\nGenerating Figure 6: Inverted repeats enrichment...")
    generate_fig6_inverted_repeats(results_dir, out_dir, format)

    # Note: Figure 7 (structured chains) is generated by analyze_structured_chains.jl

    println("\n=== Figures Complete ===")
end

main()
