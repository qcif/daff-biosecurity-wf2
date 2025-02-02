#!/usr/bin/env bash

singularity exec \
    docker://neoformit/daff-taxonomic-assignment \
    python taxonkit.py \
    test-data/taxids.csv
