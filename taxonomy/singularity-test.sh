#!/usr/bin/env bash

singularity exec \
    docker://neoformit/daff-taxonomic-assignment \
    -B .:/working \
    -B /home/ubuntu/.taxonkit:/root/.taxonkit \
    python /app/taxonomy/ncbi_taxonomy.py \
    test-data/blast-output.txt
