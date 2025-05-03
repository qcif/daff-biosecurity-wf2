#!/usr/bin/env bash

# Run a series of end-to-end tests over all workflow components, on a series of
# test cases that cover different scenarios.

set -e

TAXDUMP_DIR="$HOME/.taxonkit"
PYTHON=$(realpath --no-symlinks $(dirname "$0")/../../venv/bin/python)
SCRIPTS_ROOT=$(realpath $(dirname "$0")/../../scripts)
TEST_CASE_ROOT=$(realpath $(dirname "$0")/../test-data/integration)
echo "Running scripts with Python interpreter: $PYTHON"
echo "Reading test cases from $TEST_CASE_ROOT"
echo "Running scripts from $SCRIPTS_ROOT"

WDIR_ROOT=$(realpath ./working)

GREEN='\033[32m'
RESET='\033[0m'

for test_case in $TEST_CASE_ROOT/*; do
    if [ -d "$test_case" ]; then
        echo ""
        echo "Running test case: $test_case"

        TEST_CASE_NAME=$(basename "$test_case")
        WDIR="$WDIR_ROOT/$TEST_CASE_NAME"
        echo "Set working directory: $WDIR"
        mkdir -p "$WDIR"
        echo "Clearing working directory..."
        rm -rf "$WDIR/*"
        echo "Copying test case files to working directory..."
        cp -r "$test_case"/* "$WDIR"

        echo "Running P0 validation..."
        $PYTHON "$SCRIPTS_ROOT/p0_validation.py" \
            --metadata_csv "$WDIR/metadata.csv" \
            --query_fasta "$WDIR/query.fasta" \
            --taxdb_dir "$TAXDUMP_DIR"
        printf "${GREEN}PASS${RESET}\n"

        echo "Running P1 BLAST parse..."
        $PYTHON "$SCRIPTS_ROOT/p1_parse_blast.py" \
            $WDIR//blast_result.xml \
            --output_dir "$WDIR"
        printf "${GREEN}PASS${RESET}\n"

        echo "Running P2 extract taxonomy..."
        $PYTHON "$SCRIPTS_ROOT/p2_extract_taxonomy.py" \
            "$WDIR/taxids.csv" \
            --output_dir "$WDIR"
        printf "${GREEN}PASS${RESET}\n"

        echo "Running P3 assign taxonomy..."
        $PYTHON "$SCRIPTS_ROOT/p3_assign_taxonomy.py" \
            $WDIR/query_001* \
            --output_dir "$WDIR"
        printf "${GREEN}PASS${RESET}\n"

        echo "Running P4 source diversity..."
        $PYTHON "$SCRIPTS_ROOT/p4_source_diversity.py" \
            $WDIR/query_001* \
            --output_dir "$WDIR"
        printf "${GREEN}PASS${RESET}\n"

        echo "Running P5 database coverage..."
        $PYTHON "$SCRIPTS_ROOT/p5_db_coverage.py" \
            $WDIR/query_001* \
            --output_dir "$WDIR"
        printf "${GREEN}PASS${RESET}\n"

        echo "Copying newick tree to working directory..."
        cp $WDIR/*.nwk $WDIR/query_001*

        echo "Running P6 report generation..."
        $PYTHON "$SCRIPTS_ROOT/p6_report.py" \
            $WDIR/query_001* \
            --output_dir "$WDIR"
        printf "${GREEN}PASS${RESET}\n"
    else
        echo "Skipping non-directory in test case root: $test_case"
    fi

done

echo ""
echo -e "${GREEN}All tests completed${RESET}\n"
echo ""
