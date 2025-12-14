# Darwin Operator Genomics

Computational framework for analyzing genomic sequences through operator algebras, supporting the preprint *"The Quaternionic Syntax of Existence"*.

## Quick Start

```bash
# Clone and enter
git clone https://github.com/Chiuratto-AI/darwin-operator-genomics
cd darwin-operator-genomics

# Full reproduction (setup + download + analyze + figures)
make reproduce
```

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
make analyze

# 4. Generate publication figures
make figures

# 5. Run tests
make test
```

## Expected Runtime

| Step | Time (200 genomes) | Disk |
|------|-------------------|------|
| setup | 2-5 min | 500 MB |
| fetch_data | 10-30 min | ~8 GB |
| analyze | 5-15 min | 100 MB |
| figures | 1-2 min | 50 MB |
| **Total** | **20-50 min** | **~10 GB** |

## Output Artifacts

After `make reproduce`:

```
results/
├── figures/
│   ├── fig1_operator_map.pdf
│   ├── fig2_symmetry_stats.pdf
│   └── fig3_quaternion_vs_baselines.pdf
├── tables/
│   └── table1_dataset_summary.csv
└── text/
    ├── results_r1_r2.md
    └── results_r3.md
```

## Scientific Background

This package implements:

1. **Dihedral group actions** on circular bacterial genomes:
   - S: cyclic shift
   - R: reverse
   - K: complement
   - RC: reverse complement

2. **Symmetry analysis**: orbit sizes and fixed points under dihedral action

3. **Quaternionic compression hypothesis**: testing whether unit quaternion representations improve next-operator prediction vs. Markov baselines

## Testing

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Verifies operator identities:
- R∘R = id, K∘K = id, S^n = id
- Dihedral relation: R∘S = S⁻¹∘R
- Relocation identities for mutations
- Non-invertibility of indels

## Citation

```bibtex
@software{agourakis2024darwin,
  author = {Agourakis, Demetrios Chiuratto},
  title = {Darwin Operator Genomics},
  year = {2024},
  url = {https://github.com/Chiuratto-AI/darwin-operator-genomics}
}
```

## Related

- [Demetrios](https://github.com/Chiuratto-AI/demetrios) - Parent language project

## License

MIT License - see [LICENSE](LICENSE)
