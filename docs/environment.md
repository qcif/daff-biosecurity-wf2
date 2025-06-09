# Environment variables

Many aspects of the analysis can be controlled through environment variables,
which are read in the [config.py](https://github.com/qcif/daff-biosecurity-wf2/blob/main/scripts/src/utils/config.py) module. Most of these are optional, but some
scripts require some variables to be present (this is documented in the main
[README.md](../README.md)).

## Credentials

- `USER_EMAIL`: used to authenticate with NCBI Entrez API if `NCBI_API_KEY` not provided. Also used to rate limit requests from different users on the same system.
- `NCBI_API_KEY`: used to authenticate with NCBI Entrez API for an increased rate limit on requests.

## Reference data

- `TAXONKIT_DATA` (default `$HOME/.taxonkit`): the directory where NCBI's taxdump files can be found.
- `ALLOWED_LOCI_FILE` (default `./loci.json`): A JSON file which describes permitted loci (barcoding regions) and their synonyms. The default [loci.json](../loci.json) could be overridden with a file in the same format in future. The `ambiguous_synonyms` listed in this file are locus synonyms which may be ambiguous in a GenBank query, for example the term "COI" may match many records that have nothing to do with the Cytochrome Oxidase gene.

## User input files

- `INPUT_FASTA_FILEPATH`: the query FASTA file containing the user's sample DNA sequences (access to this is required throughout the modules)
- `INPUT_METADATA_CSV_FILEPATH`: the metadata CSV file containing the user's sample metadata (access to this is required throughout the modules)

## Working directories

- `OUTPUT_DIR`: where output data should be written to for the Nextflow run
- `QUERY_DIR`: the directory containing output data for the querying currently being analysed.

## Parameters

- `BLAST_MAX_TARGET_SEQS` (default `2000`): the maximum number of hits that was collected for each query sequence in the BLAST search. This is not used for analysis but is rendered in the report.
- `DB_COVERAGE_TOI_LIMIT` (default `10`): the maximum number of TOIs that will be analysed by P5 (database coverage).
- `GBIF_LIMIT_RECORDS` (default `500`): the maximum number of records per-request to the GBIF API. More records than this will be fetched in batches.
- `GBIF_MAX_OCCURRENCE_RECORDS` (default `5000`): the maximum number of GBIF records that will be fetched for plotting the occurrence distribution map.
- `GBIF_ACCEPTED_STATUS` (default `accepted,doubtful`): only GBIF records with these statuses will be retained when fetching related species (comma-separated)
- `BLAST_DATABASE_NAME` (default `NCBI Core Nt`): for showing in the report
- `FACILITY_NAME`: this will be shown in the report
- `ANALYST_NAME`: this will be shown in the report
- `REPORT_DEBUG`: if set to `1` this replaces the timestamp in the report file name with `DEBUG` so that it can be easily reloaded in the browser after re-rendering.

## Filtering criteria

- `MIN_NT`: minimum alignment length for a BLAST hit to be considered for candidate screening (nucleotides).
- `MIN_Q_COVERAGE` (default `300`): minimum query coverage for a BLAST hit to be considered for candidate screening (nucleotides).
- `MIN_IDENTITY` (default `0.935`): minimum identity for a BLAST hit to be considered for candidate screening (decimal).
- `MIN_IDENTITY_STRICT` (default ` 0.985`): minimum hit identity to be considered a STRONG candidate (decimal).
- `MAX_CANDIDATES_FOR_ANALYSIS` (default `3`): the maximum number of candidate species that will proceed to further analysis (P4/5). When this threshold is reached, a boxplot showing identity distributions is shown in the report "Candidates" section.
- `MIN_SOURCE_COUNT` (default `5`): Minimum number of independent publications required for a candidate species to receive Flag 4A.
- `DB_COV_MIN_A` (default `5`): Minimum number of GenBank records to receive Flag 5.1A
- `DB_COV_MIN_B` (default `1`): Minimum number of GenBank records to receive Flag 5.1B
- `DB_COV_RELATED_MIN_A` (default `90`): Minimum percent species coverage of GenBank records to receive Flag 5.2A
- `DB_COV_RELATED_MIN_B` (default `10`): Minimum percent species coverage of GenBank records to receive Flag 5.2B
- `DB_COV_COUNTRY_MISSING_A` (default `1`): Minimum number of species WITHOUT GenBank records to receive Flag 5.2B

## Output file names

These variables allow customization of output file names, if the default file name is not acceptable for some reason.

- `TIMESTAMP_FILENAME` (default `timestamp.txt`)
- `ACCESSIONS_FILENAME` (default `accessions.txt`)
- `TAXONOMY_FILENAME` (default `taxonomy.csv`)
- `QUERY_TITLE_FILENAME` (default `query_title.txt`)
- `HITS_JSON_FILENAME` (default `hits.json`)
- `HITS_FASTA_FILENAME` (default `hits.fasta`)
- `TAXONOMY_ID_CSV_FILENAME` (default `assigned_taxonomy.csv`)
- `CANDIDATES_FASTA_FILENAME` (default `candidates.fasta`)
- `CANDIDATES_CSV_FILENAME` (default `candidates.csv`)
- `CANDIDATES_JSON_FILENAME` (default `candidates.json`)
- `CANDIDATES_COUNT_FILENAME` (default `candidates_count.txt`)
- `CANDIDATES_SOURCES_JSON_FILENAME` (default `candidates_sources.json`)
- `INDEPENDENT_SOURCES_JSON_FILENAME` (default `aggregated_sources.json`)
- `TOI_DETECTED_CSV_FILENAME` (default `taxa_of_concern_detected.csv`)
- `PMI_MATCH_CSV_FILENAME` (default `preliminary_id_match.csv`)
- `BOXPLOT_IMG_FILENAME` (default `boxplot.png`)
- `TREE_NWK_FILENAME` (default `candidates.nwk`)
- `DB_COVERAGE_JSON_FILENAME` (default `db_coverage.json`)
- `BOLD_TAXON_COUNT_JSON` (default `bold_taxon_counts.json`)
- `BOLD_TAXON_COLLECTORS_JSON` (default `bold_taxon_collectors.json`)
- `BOLD_TAXONOMY_JSON` (default `bold_taxonomy.json`)