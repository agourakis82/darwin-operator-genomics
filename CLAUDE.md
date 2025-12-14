# CLAUDE.md - Repository Guidance

## Project Overview
Darwin Operator Genomics: Julia implementation of operator algebra for genomic sequence analysis, supporting the preprint "The Quaternionic Syntax of Existence".

## Key Commands
```bash
# Setup (first time)
make setup

# Download bacterial genomes (default: 200 assemblies)
make fetch_data

# Run all analyses
make analyze

# Generate figures
make figures

# Run tests
make test

# Full reproducibility (setup + fetch + analyze + figures)
make reproduce

# Clean generated outputs
make clean
```

## Directory Structure
- `src/` - Main Julia module (DarwinOperatorGenomics.jl)
- `scripts/` - CLI scripts (fetch, analyze, quaternion test)
- `data/cache/` - Downloaded genomes (gitignored)
- `results/figures/` - Generated figures
- `results/tables/` - Generated tables
- `results/text/` - Markdown snippets for paper
- `paper/` - Manuscript files
- `test/` - Unit tests

## Code Style
- Julia standard style (4-space indent)
- Use BioSequences.jl types for DNA sequences
- Document public functions with docstrings
- Operators: S (shift), R (reverse), K (complement), RC (reverse_complement)

## Testing
```bash
julia --project -e 'using Pkg; Pkg.test()'
```

## Data Manifest
All downloads logged to `data/cache/manifest.jsonl` with:
- accession, URL, timestamp, size, SHA256 checksum

## Related Projects
- [Demetrios](https://github.com/Chiuratto-AI/demetrios) - Parent language project
