# Environment variables

Many aspects of the analysis can be controlled through environment variables,
which are read in the [config.py](https://github.com/qcif/daff-biosecurity-wf2/blob/main/scripts/src/utils/config.py) module. Most of these are optional, but some
scripts require some variables to be present (this is documented in the main
[README.md](../README.md)).

<table>
    <thead>
        <tr>
            <th>Variable</th>
            <th>Default</th>
            <th>Description</th>
        </tr>
    </thead>
    <tbody>
        <tr><td colspan=3 style="background: #99999911; text-align: center;">Credentials</td></tr>
        <tr>
            <td><strong>USER_EMAIL</strong></td>
            <td></td>
            <td>Used to authenticate with NCBI Entrez API if <code>NCBI_API_KEY</code> not provided. Also used to rate limit requests from different users on the same system.</td>
        </tr>
        <tr>
            <td><strong>NCBI_API_KEY</strong></td>
            <td></td>
            <td>Used to authenticate with NCBI Entrez API for an increased rate limit on requests.</td>
        </tr>
        <tr><td colspan=3 style="background: #99999911; text-align: center;">Reference data</td></tr>
        <tr>
            <td><strong>TAXONKIT_DATA</strong></td>
            <td>$HOME/.taxonkit</td>
            <td>The directory where NCBI's taxdump files can be found.</td>
        </tr>
        <tr>
            <td><strong>ALLOWED_LOCI_FILE</strong></td>
            <td>./loci.json</td>
            <td>A JSON file which describes permitted loci (barcoding regions) and their synonyms. The default can be overridden with a file in the same format. The <code>ambiguous_synonyms</code> listed in this file are locus synonyms which may be ambiguous in a GenBank query.</td>
        </tr>
        <tr><td colspan=3 style="background: #99999911; text-align: center;">User inputs</td></tr>
        <tr>
            <td><strong>INPUT_FASTA_FILEPATH</strong></td>
            <td></td>
            <td>The query FASTA file containing the user's sample DNA sequences (required throughout the modules).</td>
        </tr>
        <tr>
            <td><strong>INPUT_METADATA_CSV_FILEPATH</strong></td>
            <td></td>
            <td>The metadata CSV file containing the user's sample metadata (required throughout the modules).</td>
        </tr>
        <tr><td colspan=3 style="background: #99999911; text-align: center;">Working directories</td></tr>
        <tr>
            <td><strong>OUTPUT_DIR</strong></td>
            <td></td>
            <td>Where output data should be written to for the Nextflow run.</td>
        </tr>
        <tr>
            <td><strong>QUERY_DIR</strong></td>
            <td></td>
            <td>The directory containing output data for the querying currently being analysed.</td>
        </tr>
        <tr><td colspan=3 style="background: #99999911; text-align: center;">Parameters</td></tr>
        <tr>
            <td><strong>BLAST_MAX_TARGET_SEQS</strong></td>
            <td>2000</td>
            <td>The maximum number of hits collected for each query sequence in the BLAST search. Not used for analysis but rendered in the report.</td>
        </tr>
        <tr>
            <td><strong>DB_COVERAGE_TOI_LIMIT</strong></td>
            <td>10</td>
            <td>The maximum number of TOIs that will be analysed by P5 (database coverage).</td>
        </tr>
        <tr>
            <td><strong>GBIF_LIMIT_RECORDS</strong></td>
            <td>500</td>
            <td>The maximum number of records per-request to the GBIF API. More records than this will be fetched in batches.</td>
        </tr>
        <tr>
            <td><strong>GBIF_MAX_OCCURRENCE_RECORDS</strong></td>
            <td>5000</td>
            <td>The maximum number of GBIF records that will be fetched for plotting the occurrence distribution map.</td>
        </tr>
        <tr>
            <td><strong>GBIF_ACCEPTED_STATUS</strong></td>
            <td>accepted,doubtful</td>
            <td>Only GBIF records with these statuses will be retained when fetching related species (comma-separated).</td>
        </tr>
        <tr>
            <td><strong>BLAST_DATABASE_NAME</strong></td>
            <td>NCBI Core Nt</td>
            <td>For showing in the report.</td>
        </tr>
        <tr>
            <td><strong>FACILITY_NAME</strong></td>
            <td></td>
            <td>This will be shown in the report.</td>
        </tr>
        <tr>
            <td><strong>ANALYST_NAME</strong></td>
            <td></td>
            <td>This will be shown in the report.</td>
        </tr>
        <tr>
            <td><strong>REPORT_DEBUG</strong></td>
            <td></td>
            <td>If set to <code>1</code> this replaces the timestamp in the report file name with <code>DEBUG</code> so that it can be easily reloaded in the browser after re-rendering.</td>
        </tr>
        <tr><td colspan=3 style="background: #99999911; text-align: center;">Analysis/filtering criteria</td></tr>
        <tr>
            <td><strong>MIN_NT</strong></td>
            <td></td>
            <td>Minimum alignment length for a BLAST hit to be considered for candidate screening (nucleotides).</td>
        </tr>
        <tr>
            <td><strong>MIN_Q_COVERAGE</strong></td>
            <td>300</td>
            <td>Minimum query coverage for a BLAST hit to be considered for candidate screening (nucleotides).</td>
        </tr>
        <tr>
            <td><strong>MIN_IDENTITY</strong></td>
            <td>0.935</td>
            <td>Minimum identity for a BLAST hit to be considered for candidate screening (decimal).</td>
        </tr>
        <tr>
            <td><strong>MIN_IDENTITY_STRICT</strong></td>
            <td>0.985</td>
            <td>Minimum hit identity to be considered a STRONG candidate (decimal).</td>
        </tr>
        <tr>
            <td><strong>MAX_CANDIDATES_FOR_ANALYSIS</strong></td>
            <td>3</td>
            <td>The maximum number of candidate species that will proceed to further analysis (P4/5). When this threshold is reached, a boxplot showing identity distributions is shown in the report "Candidates" section.</td>
        </tr>
        <tr>
            <td><strong>MIN_SOURCE_COUNT</strong></td>
            <td>5</td>
            <td>Minimum number of independent publications required for a candidate species to receive Flag 4A.</td>
        </tr>
        <tr>
            <td><strong>DB_COV_MIN_A</strong></td>
            <td>5</td>
            <td>Minimum number of GenBank records to receive Flag 5.1A.</td>
        </tr>
        <tr>
            <td><strong>DB_COV_MIN_B</strong></td>
            <td>1</td>
            <td>Minimum number of GenBank records to receive Flag 5.1B.</td>
        </tr>
        <tr>
            <td><strong>DB_COV_RELATED_MIN_A</strong></td>
            <td>90</td>
            <td>Minimum percent species coverage of GenBank records to receive Flag 5.2A.</td>
        </tr>
        <tr>
            <td><strong>DB_COV_RELATED_MIN_B</strong></td>
            <td>10</td>
            <td>Minimum percent species coverage of GenBank records to receive Flag 5.2B.</td>
        </tr>
        <tr>
            <td><strong>DB_COV_COUNTRY_MISSING_A</strong></td>
            <td>1</td>
            <td>Minimum number of species WITHOUT GenBank records to receive Flag 5.2B.</td>
        </tr>
        <tr><td colspan=3 style="background: #99999911; text-align: center;">Output file names</td></tr>
        <tr>
            <td><strong>TIMESTAMP_FILENAME</strong></td>
            <td>timestamp.txt</td>
            <td>Customize output file name if the default is not acceptable.</td>
        </tr>
        <tr>
            <td><strong>ACCESSIONS_FILENAME</strong></td>
            <td>accessions.txt</td>
            <td>Customize output file name if the default is not acceptable.</td>
        </tr>
        <tr>
            <td><strong>TAXONOMY_FILENAME</strong></td>
            <td>taxonomy.csv</td>
            <td>Customize output file name if the default is not acceptable.</td>
        </tr>
        <tr>
            <td><strong>QUERY_TITLE_FILENAME</strong></td>
            <td>query_title.txt</td>
            <td>Customize output file name if the default is not acceptable.</td>
        </tr>
        <tr>
            <td><strong>HITS_JSON_FILENAME</strong></td>
            <td>hits.json</td>
            <td>Customize output file name if the default is not acceptable.</td>
        </tr>
        <tr>
            <td><strong>HITS_FASTA_FILENAME</strong></td>
            <td>hits.fasta</td>
            <td>Customize output file name if the default is not acceptable.</td>
        </tr>
        <tr>
            <td><strong>TAXONOMY_ID_CSV_FILENAME</strong></td>
            <td>assigned_taxonomy.csv</td>
            <td>Customize output file name if the default is not acceptable.</td>
        </tr>
        <tr>
            <td><strong>CANDIDATES_FASTA_FILENAME</strong></td>
            <td>candidates.fasta</td>
            <td>Customize output file name if the default is not acceptable.</td>
        </tr>
        <tr>
            <td><strong>CANDIDATES_CSV_FILENAME</strong></td>
            <td>candidates.csv</td>
            <td>Customize output file name if the default is not acceptable.</td>
        </tr>
        <tr>
            <td><strong>CANDIDATES_JSON_FILENAME</strong></td>
            <td>candidates.json</td>
            <td>Customize output file name if the default is not acceptable.</td>
        </tr>
        <tr>
            <td><strong>CANDIDATES_COUNT_FILENAME</strong></td>
            <td>candidates_count.txt</td>
            <td>Customize output file name if the default is not acceptable.</td>
        </tr>
        <tr>
            <td><strong>CANDIDATES_SOURCES_JSON_FILENAME</strong></td>
            <td>candidates_sources.json</td>
            <td>Customize output file name if the default is not acceptable.</td>
        </tr>
        <tr>
            <td><strong>INDEPENDENT_SOURCES_JSON_FILENAME</strong></td>
            <td>aggregated_sources.json</td>
            <td>Customize output file name if the default is not acceptable.</td>
        </tr>
        <tr>
            <td><strong>TOI_DETECTED_CSV_FILENAME</strong></td>
            <td>taxa_of_concern_detected.csv</td>
            <td>Customize output file name if the default is not acceptable.</td>
        </tr>
        <tr>
            <td><strong>PMI_MATCH_CSV_FILENAME</strong></td>
            <td>preliminary_id_match.csv</td>
            <td>Customize output file name if the default is not acceptable.</td>
        </tr>
        <tr>
            <td><strong>BOXPLOT_IMG_FILENAME</strong></td>
            <td>boxplot.png</td>
            <td>Customize output file name if the default is not acceptable.</td>
        </tr>
        <tr>
            <td><strong>TREE_NWK_FILENAME</strong></td>
            <td>candidates.nwk</td>
            <td>Customize output file name if the default is not acceptable.</td>
        </tr>
        <tr>
            <td><strong>DB_COVERAGE_JSON_FILENAME</strong></td>
            <td>db_coverage.json</td>
            <td>Customize output file name if the default is not acceptable.</td>
        </tr>
        <tr>
            <td><strong>BOLD_TAXON_COUNT_JSON</strong></td>
            <td>bold_taxon_counts.json</td>
            <td>Customize output file name if the default is not acceptable.</td>
        </tr>
        <tr>
            <td><strong>BOLD_TAXON_COLLECTORS_JSON</strong></td>
            <td>bold_taxon_collectors.json</td>
            <td>Customize output file name if the default is not acceptable.</td>
        </tr>
        <tr>
            <td><strong>BOLD_TAXONOMY_JSON</strong></td>
            <td>bold_taxonomy.json</td>
            <td>Customize output file name if the default is not acceptable.</td>
        </tr>
    </tbody>
</table>
