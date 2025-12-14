#!/usr/bin/env julia
"""
One-command reproduction script.

Runs the complete analysis pipeline:
1. Downloads genomes (or uses cache)
2. Runs symmetry analysis
3. Runs quaternion test
4. Generates figures
5. Updates paper assets
"""

using Dates

function run_step(name::String, cmd::Cmd)
    println("\n" * "="^60)
    println("STEP: $name")
    println("="^60 * "\n")

    t0 = time()
    result = run(cmd; wait=true)
    elapsed = round(time() - t0, digits=1)

    println("\n✓ Completed in $(elapsed)s\n")
    return result
end

function main()
    println("""

    ╔══════════════════════════════════════════════════════════╗
    ║     Darwin Operator Genomics - Reproduction Script       ║
    ║                                                          ║
    ║  This script regenerates all results from scratch.       ║
    ╚══════════════════════════════════════════════════════════╝

    Started: $(now())

    """)

    t_total = time()

    # Parse command line arguments
    max_genomes = get(ENV, "MAX", "200")
    seed = get(ENV, "SEED", "42")

    println("Configuration:")
    println("  MAX genomes: $max_genomes")
    println("  SEED: $seed")

    # Step 1: Fetch data
    run_step("Fetch NCBI Complete Bacterial Genomes",
        `julia --project=. scripts/fetch_ncbi_complete_bacteria.jl --max $max_genomes --seed $seed --resume`)

    # Step 2: Analyze genomes
    run_step("Analyze Genome Symmetry Statistics",
        `julia --project=. scripts/analyze_genomes.jl`)

    # Step 3: Generate figures 1 & 2
    run_step("Generate Publication Figures",
        `julia --project=. scripts/generate_figures.jl`)

    # Step 4: Quaternion test
    run_step("Quaternionic Compression Hypothesis Test",
        `julia --project=. scripts/quaternion_test.jl`)

    # Step 5: Update paper
    run_step("Update Paper Assets",
        `bash -c "mkdir -p paper/figures && cp -f results/figures/*.pdf paper/figures/"`)

    # Summary
    total_time = round((time() - t_total) / 60, digits=1)

    println("""

    ╔══════════════════════════════════════════════════════════╗
    ║                   REPRODUCTION COMPLETE                  ║
    ╚══════════════════════════════════════════════════════════╝

    Total time: $(total_time) minutes

    Generated artifacts:
    - results/tables/table1_dataset_summary.csv
    - results/tables/replicon_stats.csv
    - results/tables/window_stats.csv
    - results/tables/quaternion_results.csv
    - results/figures/fig1_operator_map.pdf
    - results/figures/fig2_symmetry_stats.pdf
    - results/figures/fig3_quaternion_vs_baselines.pdf
    - results/text/results_r1_r2.md
    - results/text/results_r3.md
    - paper/figures/*.pdf

    To view results:
      cat results/text/results_r1_r2.md
      cat results/text/results_r3.md

    Finished: $(now())
    """)
end

main()
