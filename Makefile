PYTHON=python3

all: help

help:
	@echo "---------------------------HELP------------------------"
	@echo "To clean files from a previous run type make clean"
	@echo "-------------------------------------------------------"

clean:
	@echo "Cleaning up output files from test environment"
	@rm -r tests/TEST_GENE/ tests/*.log tests/percent_carriers_in_population.txt