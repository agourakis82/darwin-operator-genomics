.PHONY: setup fetch_data analyze figures paper test reproduce clean help

JULIA := julia --project=.
MAX ?= 200
SEED ?= 42

help:
	@echo "Darwin Operator Genomics - Makefile targets"
	@echo ""
	@echo "  setup       Install Julia dependencies"
	@echo "  fetch_data  Download NCBI bacterial genomes (MAX=$(MAX), SEED=$(SEED))"
	@echo "  analyze     Run genome symmetry analysis"
	@echo "  figures     Generate all publication figures"
	@echo "  paper       Update paper assets from results"
	@echo "  test        Run unit tests"
	@echo "  reproduce   Full pipeline: setup + fetch + analyze + figures + paper"
	@echo "  clean       Remove generated outputs (keeps cache)"
	@echo "  cleanall    Remove everything including cache"
	@echo ""
	@echo "Examples:"
	@echo "  make reproduce          # Full run with defaults"
	@echo "  make fetch_data MAX=50  # Smaller dataset"

setup:
	$(JULIA) -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

fetch_data:
	$(JULIA) scripts/fetch_ncbi_complete_bacteria.jl --max $(MAX) --seed $(SEED) --out data/cache/

analyze:
	$(JULIA) scripts/analyze_genomes.jl --cache data/cache/ --out results/

figures: analyze
	$(JULIA) scripts/generate_figures.jl --results results/ --out results/figures/

quaternion_test:
	$(JULIA) scripts/quaternion_test.jl --cache data/cache/ --out results/

paper: figures quaternion_test
	@mkdir -p paper/figures
	@cp -f results/figures/*.pdf paper/figures/ 2>/dev/null || true
	@cp -f results/figures/*.png paper/figures/ 2>/dev/null || true
	@echo "Paper assets updated in paper/figures/"

test:
	$(JULIA) -e 'using Pkg; Pkg.test()'

reproduce: setup fetch_data analyze figures quaternion_test paper
	@echo ""
	@echo "============================================"
	@echo "Reproduction complete!"
	@echo "============================================"
	@echo "Figures: results/figures/"
	@echo "Tables:  results/tables/"
	@echo "Text:    results/text/"
	@echo "Paper:   paper/figures/"
	@echo "============================================"

clean:
	rm -rf results/figures/* results/tables/* results/text/*
	rm -rf paper/figures/*

cleanall: clean
	rm -rf data/cache/*
