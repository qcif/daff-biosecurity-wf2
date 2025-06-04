# DAFF Biosecurity Workflows - Taxonomic Assignment

These Python modules are used for intermediate data processing as part of a
[Nextflow workflow](https://github.com/qcif/nf-daff-biosecurity-wf2/tree/main).
A user guide for each module is listed here. Examples of running
each script can be found in the [.vscode/launch.json](.vscode/launch.json) file,
which shows CLI arguments and environment variables for each script.

Throughout execution of these scripts, access to the inputs files is required.
To avoid repeated passing of these files, they can just be set as environment
variables:

```sh
INPUT_FASTA_FILEPATH="/my/input/folder/query.fasta"
INPUT_METADATA_CSV_FILEPATH="/my/input/folder/metadata.csv"
LOGGING_DEBUG=0  # 1 to enable
SKIP_ORIENTATION=0  # 1 to skip orientation of BOLD sequences (requires setup)
GBIF_MAX_OCCURRENCE_RECORDS=200  # Reduce to 200 for testing/dev to speed up p5. Default 5000.
REPORT_DEBUG=0  # 1 to omit timestamp from report filename for browser reload between changes
FACILITY_NAME="Hogwarts"  # Displayed in report
ANALYST_NAME="Harry Potter"  # Displayed in report
```

## P0 validate inputs

This script takes user inputs and validates them to ensure that the format and
content are valid for the requested analysis. If one of the inputs are found to
have erroneous content, this will raise an exception with an error message that
should be sufficient for the user to understand what's wrong with their input
data.

See [config.INPUTS](https://github.com/qcif/daff-biosecurity-wf2/blob/main/scripts/src/utils/config.py#L113)
for some parameters which are used for validation, such as permitted FASTA
sequence lengths and required metadata.csv fields.

```
$ python p0_validation.py -h

usage: p0_validation.py [-h] --metadata_csv METADATA_CSV --query_fasta
                        QUERY_FASTA --taxdb_dir TAXDB_DIR [--bold]

Validate user input.

options:
  -h, --help            show this help message and exit
  --metadata_csv METADATA_CSV
                        Path to metadata.csv input file.
  --query_fasta QUERY_FASTA
                        Path to queries.fasta input file.
  --taxdb_dir TAXDB_DIR
                        Path to queries.fasta input file.
  --bold                Validate inputs for a BOLD analysis (accept blank
                        locus field).
```


## P1 BLAST parser

This script parses a BLAST XML output into JSON format that is more easily
readable by Python and other programs. It also extracts FASTA sequences for
BLAST hit subjects, and writes a list of hit accessions to be used as the input
for the next workflow step.

```
$ python scripts/p1_parse_blast.py -h

Parse BLAST XML output file.

positional arguments:
  blast_xml_path        Path to the BLAST XML file to parse.

options:
  -h, --help            show this help message and exit
  --input-db INPUT_DB   Database path to use for retrieving taxon ID.
  --output_dir OUTPUT_DIR
                        Directory to save parsed output files (JSON and FASTA).
```

Output is in per-query directories corresponding to the sequence index in
the query FASTA file:

```
output/
├── accessions.txt  # input file for blastdbcmd
├── query_001
│   ├── hits.fasta
│   ├── hits.json
│   └── query_title.txt
├── query_002
│   ├── hits.fasta
│   ├── hits.json
│   └── query_title.txt
├── ...
...
```

## BLASTDBCMD

This is not a Python module but is a required intermediate step between modules
1 and 2 implemented in the parent Nextflow workflow. The `blastdbcmd` tool
should be used to extract taxids for each accession as follows:

```sh
blastdbcmd -entry_batch output/accessions.txt -db </path/to/blastdb> -outfmt "%a,%T" > output/taxids.csv
```

NOTE: Some weird behaviour observed by this tool - it seems to extract more
accessions than it is provided. I passed a file of 1830 taxids and found that
2169 were written. There were no duplicates, just extra accessions that weren't
in the input list! For the purpose of this pipeline it's a no-op, but worth
noting the unexpected behaviour (bug).


## P2 NCBI Taxonomy extractor

This script is used for fetching taxonomy information for a list of taxids. The
input should be provided in CSV format with columns (accession,taxid). The
output is a CSV file with columns (accession,taxid,superkingdom,kingdom,phylum,
class,order,family,genus,species). This script requires access to both a
`taxonkit` executable in PATH and NCBI taxdump directory. A custom directory
for the latter can be specified with the `--taxdb` param.

```
$ python scripts/p2_extract_taxonomy.py -h
usage: ncbi_taxonomy.py [-h] [--output OUTPUT_CSV] [--taxdb TAXDB_PATH] taxids_csv

Extract taxids and taxonomic information from NCBI databases. This requires access to the NCBI taxdump files via a CLI argument.

positional arguments:
  taxids_csv           CSV file with columns (accession,taxid) to extract taxonomy information for.

options:
  -h, --help           show this help message and exit
  --output OUTPUT_CSV  CSV file where taxonomy data will be written. Defaults to taxonomy.csv
  --taxdb TAXDB_PATH   Path to directory containing NCBI taxdump files for taxonkit. Defaults to /home/ubuntu/.taxonkit
  ```


## P3 Evaluate taxonomy

This script evaluates BLAST/BOLD hits enumerating the number of hits that fall
within the candidate match criteria, and then enumerating the number of species
represented by those hits. It generates several reportable outcomes:

- Attempt to assign a taxonomic identity for the sample (Flag 1)
- Check for Taxa of Interest matching the sample (Flag 2)
- Check whether the Preliminary ID matches the taxonomic identity (Flag 7; only if Flag 1A)
- If more than 3 candidate species, produces a boxplot showing distribution of hit identity % for each species

```
$ python p3_assign_taxonomy.py -h
usage: p3_assign_taxonomy.py [-h] [--output_dir OUTPUT_DIR] [--bold] query_dir

Run logic for pipeline phase 1-2. - Attempt species ID from BLAST results.json (flag 1) - Detect Taxa of Interest (flag 2) Taxa of Interest output has the following CSV
fields: Taxon of interest: The provided TOI that matched a candidate species Match rank: The taxonomic rank of the match (TOI rank may be above species) Match taxon: The taxon
that matched the TOI Match species: The species of the candidate match Match accession: The NCBI accession of the candidate match

positional arguments:
  query_dir             Path to query output directory

options:
  -h, --help            show this help message and exit
  --output_dir OUTPUT_DIR
                        Path to output directory. Defaults to output.
  --bold                Outputs are from BOLD query.
```


## P4 Analysis of reference sequence publications

This script fetches publications from GenBank metadata for each BLAST/BOLD hit
for each candidate species. The publications are then grouped into "independent
publications" depending on whether they are from the same authors. For BOLD hits
which lack a GenBank record, the BOLD `collector` field is used instead. This
analysis results in an aggregated_sources.json file which enumerates the
independent publication sources for each candidate species by hit, and lists the
publication author, journal and title for each publication.

```
$ python p4_source_diversity.py -h

usage: p4_source_diversity.py [-h] [--output_dir OUTPUT_DIR] query_dir

Analyze the diversity of reference sequence sources oer target. A source is
defined as a publication or set of authors that are linked to the genbank
record for that sequence. If there are no references, no sources are returned
and the sequence is classified as "anonymous". Many anonymous records are from
automated genome annotation projects, often carried out by NCBI themselves.
These records are flagged so that the user can be aware of the potential
reduced credibility of these annotation.

positional arguments:
  query_dir             Path to query output directory

options:
  -h, --help            show this help message and exit
  --output_dir OUTPUT_DIR
                        Path to output directory. Defaults to output.
```


## P5 Analysis of database coverage

This is the most complex analysis, because it forks multiple times into three
different analyses (5.1, 5.2, 5.3) that are run against each target taxon in
the candidates, PMI and TOIs. This results in a minimum of 3x1=3 analyses (no
candidates, no TOI provided) and beyond 3x6=18 analyses (3 candidates, 2 TOIs
provided). If more than 3 candidate species are present, they will be dropped
from the targets and the analysis will only be performed for the PMI and TOIs.
This analysis also includes the generation of geographic occurrence maps for
each target taxon.

The code for these analyses is fairly abstract in order to accomodate the range
of cases described above, and is heavily threaded to complete the analysis in a
reasonable time. The [Throttle](#throttling-api-requests)
is critical here to avoid exceeding API rate limits, since the analysis involves
sending LOTS of API requests (sometimes many hundreds per-sample). The entrypoint
for this analysis is
[assess.py](https://github.com/qcif/daff-biosecurity-wf2/blob/main/scripts/src/coverage/assess.py).

1. **Setup** (see [targets.py](https://github.com/qcif/daff-biosecurity-wf2/blob/main/scripts/src/coverage/targets.py))
  1. A list of target taxa is generated from candidate species, PMIs and TOIs
  1. TaxIDs are extracted for each target taxon
  1. GBIF records are extracted for each target taxon
1. **Occurrence maps** are drawn (see [maps.py](https://github.com/qcif/daff-biosecurity-wf2/blob/main/scripts/src/gbif/maps.py))
1. **Generate tasks**: a list of analysis tasks (targets x 3 analyses) is generated for threading
1. **Thread tasks** - for each target taxon:
  1. **5.1** - DB coverage of target taxon. How many records are in the
    reference database for the target taxon?
  1. **5.2** - DB coverage of species in target genus (only applies to targets
    at rank genus or species).
  1. **5.3** - DB coverage of species in target genus, limited to the sample
    country of origin (declared in metadata.csv input)
1. Collect results from threads and write to `db_coverage.json`


## Throttling API requests

The P5 script results in many threads across multiple independent processes
(one per-sample, as orchestrated by Nextflow) sending LOTS of API requests.
Without some kind of throttling mechanism, we would quickly get blocked by NCBI,
GBIF, BOLD, or other APIs that we access in the workflow. For this reason, ALL
requests to external APIs need to go through the Throttle. The example below
shows the Throttle being used via its `with_retry()` method to retry the
request several times before raising an exception.

```py
kwargs = {
    'species': 'Homo sapiens',
    'country': 'CA',
}
throttle = Throttle(ENDPOINTS.GBIF_FAST)
res = throttle.with_retry(
    pygbif.occurrences.search,
    kwargs=kwargs,
)
```

The throttle works by writing timestamps to a SQLite3 table which is stored
in the user's temp files and shared across instances. When a thread is being
throttled, the timestamps in the database are being checked until they indicate
that less than 10 requests have been sent in the last second. At that point,
whichever thread is "lucky enough" to acquire the database lock first is able to
write its timestamp to the table and then the throttle for that thread is
released. The database table depends on the `ENDPOINT` that the throttle was
created with:

https://github.com/qcif/daff-biosecurity-wf2/blob/cfb1a93917cb91dc81bccf8b4956a7a081a86a29/scripts/src/utils/throttle.py#L14-L30
