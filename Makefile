PYTHON=python3

all: help

help:
	@echo "---------------------------HELP------------------------"
	@echo "To clean files from a previous run type make clean"
	@echo "-------------------------------------------------------"

clean:
	@echo "Cleaning up output files from test environment"
	@rm tests/*.txt tests/*.log