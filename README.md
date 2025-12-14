# Darwin Operator Genomics

<!-- Zenodo badge placeholder - update after first release -->
<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX) -->

Computational framework for analyzing genomic sequences through operator algebras, supporting the preprint *"The Quaternionic Syntax of Existence: Operator Algebras for Genomic Evolution"*.

## Quick Start

```bash
# Clone and enter
git clone https://github.com/agourakis82/darwin-operator-genomics
cd darwin-operator-genomics

# Full reproduction (setup + download + analyze + figures)
make reproduce
```

## Features

- **Exact symmetry analysis**: Orbit sizes and fixed points under dihedral action
- **Approximate symmetry**: d_min/L metric with GC-shuffled null baseline
- **Binary dihedral representation**: Dicyclic group embedding in unit quaternions
- **Quaternion compression test**: Comparison against Markov baselines
- **Full reproducibility**: Single-command pipeline with deterministic caching

## Requirements

- Julia 1.9+
- ~10 GB disk space for genome cache (default 200 assemblies)
- Internet connection for NCBI downloads

## Manual Steps

```bash
# 1. Install Julia dependencies
make setup

# 2. Download bacterial genomes (configurable)
make fetch_data           # Default: 200 complete genomes
make fetch_data MAX=50    # Smaller dataset

# 3. Run analyses
make analyze              # Exact symmetry (orbit sizes, fixed points)
make analyze_approx       # Approximate symmetry (d_min/L metric)

# 4. Run quaternion experiments
make quaternion_test      # Quaternion vs Markov + binary dihedral tests

# 5. Generate publication figures
make figures

# 6. Run tests
make test
```

## Expected Runtime

| Step | Time (200 genomes) | Disk |
|------|-------------------|------|
| setup | 2-5 min | 500 MB |
| fetch_data | 10-30 min | ~8 GB |
| analyze | 5-15 min | 100 MB |
| analyze_approx | 10-20 min | 50 MB |
| quaternion_test | 2-5 min | 10 MB |
| figures | 1-2 min | 50 MB |
| **Total** | **30-75 min** | **~10 GB** |

## Output Artifacts

After `make reproduce`:

```
results/
├── figures/
│   ├── fig1_operator_map.pdf
│   ├── fig2_symmetry_stats.pdf      # Exact + approx symmetry panels
│   └── fig3_quaternion_vs_baselines.pdf  # Quaternion + binary dihedral
├── tables/
│   ├── table1_dataset_summary.csv
│   ├── window_stats.csv
│   ├── replicon_stats.csv
│   ├── approx_symmetry_stats.csv
│   ├── approx_symmetry_summary.csv
│   ├── quaternion_results.csv
│   ├── dicyclic_cover_verification.csv
│   ├── dicyclic_chain_comparison.csv
│   └── dicyclic_trajectory_stats.csv
└── text/
    ├── results_r1_r2.md
    ├── results_approx_symmetry.md
    └── results_r3.md
```

## Scientific Background

This package implements:

### 1. Dihedral Group Actions (D_n)
Operators on circular bacterial genomes:
- **S**: cyclic shift (rotation)
- **R**: reverse
- **K**: complement (A↔T, C↔G)
- **RC**: reverse complement

### 2. Exact Symmetry Analysis
- Orbit sizes under dihedral action
- Fixed point detection (palindromes, RC-palindromes)
- Canonical representation (lexicographically minimal)

### 3. Approximate Dihedral Self-Similarity
- **d_min/L**: Minimum normalized Hamming distance to any dihedral transform
- GC-preserving shuffled baseline for null hypothesis
- KS test and Cohen's d effect size

### 4. Binary Dihedral (Dicyclic) Representation
- Lift D_n to Dic_n (dicyclic group of order 4n)
- Embedding in unit quaternions SU(2)
- Verified 2-to-1 double cover property
- Reference: Conway & Smith, "On Quaternions and Octonions" (2003)

### 5. Quaternion Compression Hypothesis
- Test quaternion representations for next-operator prediction
- Comparison against Markov(1) and Markov(2) baselines
- Negative result preserved for scientific rigor

## Testing

```bash
make test
# or
julia --project=. test/runtests.jl
```

Verifies:
- Operator identities: R∘R = id, K∘K = id, S^n = id
- Dihedral relation: R∘S = S⁻¹∘R
- Relocation identities for mutations
- Non-invertibility of indels
- Binary dihedral group structure and double cover property

## Zenodo Archive

To create a Zenodo-ready archive:

```bash
make zenodo_bundle
# Creates: dist/darwin-operator-genomics-YYYYMMDD.tar.gz
```

## Citation

```bibtex
@software{agourakis2024darwin,
  author = {Agourakis, Demetrios Chiuratto},
  title = {Darwin Operator Genomics: Operator Algebras for Genomic Evolution},
  year = {2024},
  url = {https://github.com/agourakis82/darwin-operator-genomics}
}
```

## References

1. Conway, J.H. & Smith, D.A. (2003). On Quaternions and Octonions. A K Peters/CRC Press.
2. O'Leary, N.A. et al. (2016). Reference sequence (RefSeq) database at NCBI. Nucleic Acids Res 44(D1):D733-D745.
3. Weyl, H. (1950). The Theory of Groups and Quantum Mechanics. Dover.

## Related

- [Demetrios](https://github.com/Chiuratto-AI/demetrios) - Parent language project

## License

MIT License - see [LICENSE](LICENSE)
