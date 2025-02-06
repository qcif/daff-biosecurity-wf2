# DAFF Biosecurity Workflows - Taxonomic Assignment

These Python modules are used for intermediate data processing as part of a
NextFlow workflow. A user guide for each module is listed here.

## BLAST parser

This script parses a BLAST XML output into JSON format that is more easily
readable by Python and other programs.

```sh
$ python blast/parse.py --help

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
├── query_1
│   ├── hits.fasta
│   ├── hits.json
│   └── query_title.txt
├── query_2
│   ├── hits.fasta
│   ├── hits.json
│   └── query_title.txt
├── query_N
...
```

## BLASTDBCMD

This is not a Python module but is a required intermediate step between modules
1 and 2. The `blastdbcmd` tool should be used to extract taxids for each
accession as follows:

```sh
blastdbcmd -entry_batch output/accessions.txt -db </path/to/blastdb> -outfmt "%a,%T" > output/taxids.csv
```

NOTE: Some weird behaviour observed by this tool - it seems to extract more
accessions than it is provided. I passed a file of 1830 taxids and found that
2169 were written. There were no duplicates, just extra accessions that weren't
in the input list! For the purpose of this pipeline it's a no-op, but worth
noting the unexpected behaviour (bug).


## NCBI Taxonomy extractor

This script is used for fetching taxonomy information for a list of taxids. The
input should be provided in CSV format with columns (accession,taxid). The
output is a CSV file with columns (accession,taxid,superkingdom,kingdom,phylum,
class,order,family,genus,species). This script requires access to both a
`taxonkit` executable in PATH and NCBI taxdump directory. A custom directory
for the latter can be specified with the `--taxdb` param.

```sh
$ python ncbi_taxonomy.py -h
usage: ncbi_taxonomy.py [-h] [--output OUTPUT_CSV] [--taxdb TAXDB_PATH] taxids_csv

Extract taxids and taxonomic information from NCBI databases. This requires access to the NCBI taxdump files via a CLI argument.

positional arguments:
  taxids_csv           CSV file with columns (accession,taxid) to extract taxonomy information for.

options:
  -h, --help           show this help message and exit
  --output OUTPUT_CSV  CSV file where taxonomy data will be written. Defaults to taxonomy.csv
  --taxdb TAXDB_PATH   Path to directory containing NCBI taxdump files for taxonkit. Defaults to /home/ubuntu/.taxonkit
  ```
