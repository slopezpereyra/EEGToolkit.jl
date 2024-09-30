# Makefile for testing
JULIA = julia
TEST_FILES = $(wildcard test/*.jl)

PHONY: all jtest clean

all:
	@echo "Building project (if necessary)..."

jtest:
	@echo "Running tests..."
	@for testfile in $(TEST_FILES); do \
		echo "Running $$testfile..."; \
		$(JULIA) $$testfile || exit 1; \
	done
	@echo "All tests passed."

# Clean target (optional)
clean:
	@echo "Cleaning up..."
