#!/usr/bin/env bash

singularity exec \
    docker://neoformit/daff-taxonomic-assignment \
    python parse.py \
    test-data/output.xml
