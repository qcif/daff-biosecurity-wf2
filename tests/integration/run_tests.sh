#!/usr/bin/env bash

# Run all integration tests from a terminal

if [[ ! -d "tests/integration" ]]; then
    echo "Integration tests directory not found. Please run from the project root."
    exit 1
fi

export PYTHONPATH="$PWD/scripts"
export TAXONKIT_DATA="$HOME/.taxonkit"
export KEEP_OUTPUTS=1   # Always retain outputs for debugging
export LOGGING_DEBUG=1
# export RUN_TEST_CASE=A  # Optional, run a specific test case

python -m unittest discover -f -v \
    -s tests/integration \
    -p test_integration.py
