#!/usr/bin/env bash

singularity exec \
    docker://neoformit/daff-taxonomic-assignment \
    -B .:/working \
    -B /home/ubuntu/.taxonkit:/root/.taxonkit \
    python /app/blast/parse.py \
    test-data/output.xml
