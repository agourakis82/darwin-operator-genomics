.PHONY: setup fetch_data analyze analyze_approx analyze_kmer analyze_gc_skew analyze_ir analyze_structured figures paper test reproduce clean cleanall help zenodo_bundle

JULIA := julia --project=.
MAX ?= 200
SEED ?= 42
WINDOWS ?= 100,250,500,1000
SAMPLES_PER_SIZE ?= 50
KMAX ?= 10
N_IR_SAMPLES ?= 5
N_STRUCT_CHAINS ?= 500

help:
	@echo "Darwin Operator Genomics - Makefile targets"
	@echo ""
	@echo "  setup         Install Julia dependencies"
	@echo "  fetch_data    Download NCBI bacterial genomes (MAX=$(MAX), SEED=$(SEED))"
	@echo "  analyze       Run genome symmetry analysis (exact invariance)"
	@echo "  analyze_approx Run approximate symmetry analysis (d_min/L)"
	@echo "  analyze_kmer  Run k-mer inversion symmetry analysis"
	@echo "  analyze_gc_skew Run GC skew and replichore analysis"
	@echo "  analyze_ir    Run inverted repeats enrichment analysis"
	@echo "  analyze_structured Run structured operator chains analysis"
	@echo "  figures       Generate all publication figures"
	@echo "  quaternion_test  Run quaternion/binary dihedral experiments"
	@echo "  paper         Update paper assets from results"
	@echo "  test          Run unit tests"
	@echo "  reproduce     Full pipeline: setup + fetch + all analyses + paper"
	@echo "  clean         Remove generated outputs (keeps cache)"
	@echo "  cleanall      Remove everything including cache"
	@echo "  zenodo_bundle Create Zenodo archive bundle"
	@echo ""
	@echo "Examples:"
	@echo "  make reproduce           # Full run with defaults"
	@echo "  make fetch_data MAX=50   # Smaller dataset"
	@echo "  make analyze_approx WINDOWS=100,500 SAMPLES_PER_SIZE=100"

setup:
	$(JULIA) -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

fetch_data:
	$(JULIA) scripts/fetch_ncbi_complete_bacteria.jl --max $(MAX) --seed $(SEED) --out data/cache/

analyze:
	$(JULIA) scripts/analyze_genomes.jl --cache data/cache/ --out results/

analyze_approx: analyze
	$(JULIA) scripts/analyze_approx_symmetry.jl --cache data/cache/ --out results/ --windows $(WINDOWS) --samples-per-size $(SAMPLES_PER_SIZE) --seed $(SEED)

analyze_kmer: analyze
	$(JULIA) scripts/analyze_kmer_symmetry.jl --cache data/cache/ --out results/ --kmax $(KMAX) --seed $(SEED)

analyze_gc_skew: analyze
	$(JULIA) scripts/analyze_gc_skew.jl --cache data/cache/ --out results/ --kmax $(KMAX)

analyze_ir: analyze
	$(JULIA) scripts/analyze_inverted_repeats.jl --cache data/cache/ --out results/ --n-samples $(N_IR_SAMPLES) --seed $(SEED)

analyze_structured:
	$(JULIA) scripts/analyze_structured_chains.jl --out results/ --n-chains $(N_STRUCT_CHAINS) --seed $(SEED)

figures: analyze_approx analyze_kmer analyze_gc_skew analyze_ir analyze_structured
	$(JULIA) scripts/generate_figures.jl --results results/ --out results/figures/

quaternion_test:
	$(JULIA) scripts/quaternion_test.jl --cache data/cache/ --out results/

paper: figures quaternion_test
	@mkdir -p paper/figures
	@cp -f results/figures/*.pdf paper/figures/ 2>/dev/null || true
	@cp -f results/figures/*.png paper/figures/ 2>/dev/null || true
	@echo "Paper assets updated in paper/figures/"

test:
	$(JULIA) test/runtests.jl

reproduce: setup fetch_data analyze analyze_approx analyze_kmer analyze_gc_skew analyze_ir analyze_structured figures quaternion_test paper
	@echo ""
	@echo "============================================"
	@echo "Reproduction complete!"
	@echo "============================================"
	@echo "Figures: results/figures/"
	@echo "Tables:  results/tables/"
	@echo "Text:    results/text/"
	@echo "Paper:   paper/figures/"
	@echo ""
	@echo "New vNext analyses:"
	@echo "  - K-mer inversion symmetry (results_kmer_symmetry.md)"
	@echo "  - GC skew & replichore (gc_skew_ori_ter.csv)"
	@echo "  - Inverted repeats (ir_enrichment_summary.csv)"
	@echo "  - Structured chains (structured_chains_summary.csv)"
	@echo "============================================"

zenodo_bundle:
	@echo "Creating Zenodo bundle..."
	@mkdir -p dist
	@tar --exclude='data/cache/*.gz' --exclude='.git' --exclude='dist' \
		-czf dist/darwin-operator-genomics-$$(date +%Y%m%d).tar.gz \
		.
	@echo "Bundle created: dist/darwin-operator-genomics-$$(date +%Y%m%d).tar.gz"
	@ls -lh dist/*.tar.gz

clean:
	rm -rf results/figures/* results/tables/* results/text/*
	rm -rf paper/figures/*

cleanall: clean
	rm -rf data/cache/*
	rm -rf dist/*
