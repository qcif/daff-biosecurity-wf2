"""Extract taxids and taxonomic information from NCBI databases.

This requires access to the NCBI taxdump files via a CLI argument.

"""

import argparse
import csv
import subprocess
import sys
import tempfile
from pathlib import Path

DEFAULT_OUTPUT_CSV = 'taxonomy.csv'
TAXONKIT_DATA = '/home/ubuntu/.taxonkit'
TAXONOMIC_RANKS = [
    "superkingdom",
    'kingdom',
    'phylum',
    'class',
    'order',
    'family',
    'genus',
    'species',
]


def main():
    args = _parse_args()
    with args.taxids_csv.open() as taxids_file:
        accession_taxids = {
            row['accession']: row['taxid']
            for row in csv.DictReader(taxids_file)
        }
    taxids = sorted(set(accession_taxids.values()))
    taxonomies = extract_taxonomies(taxids, taxdb=args.taxdb_path)
    _write_csv(taxonomies, accession_taxids, args.output_csv)


def _parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'taxids_csv',
        type=Path,
        help='CSV file with columns (accession,taxid) to extract taxonomy'
             ' information for.',
    )
    parser.add_argument(
        '--output',
        dest='output_csv',
        type=Path,
        help='CSV file where taxonomy data will be written. Defaults to'
             f' {DEFAULT_OUTPUT_CSV}',
        default=DEFAULT_OUTPUT_CSV,
    )
    parser.add_argument(
        '--taxdb',
        dest='taxdb_path',
        type=Path,
        help='Path to directory containing NCBI taxdump files for taxonkit.'
             f' Defaults to {TAXONKIT_DATA}',
        default=Path(TAXONKIT_DATA),
    )
    return parser.parse_args()


def _write_csv(taxonomies, accession_taxids, output_csv):
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open('w') as output_file:
        writer = csv.DictWriter(
            output_file,
            fieldnames=['accession', 'taxid'] + TAXONOMIC_RANKS,
        )
        writer.writeheader()
        rows = [
            {
                'accession': accession,
                'taxid': taxid,
                **taxonomies[taxid]
            }
            for accession, taxid in accession_taxids.items()
            if taxid in taxonomies
        ]
        writer.writerows(rows)


def extract_taxonomies(
    taxids: list[str],
    taxdb: Path = Path(TAXONKIT_DATA),
) -> dict[str, dict[str, str]]:
    """Use taxonkit lineage to extract taxonomic data for given taxids."""
    with tempfile.NamedTemporaryFile(mode='w+') as temp_file:
        temp_file.write("\n".join(taxids))
        temp_file.flush()
        temp_file_name = temp_file.name
        result = subprocess.run(
            [
                'taxonkit',
                'lineage',
                '-R',
                '-c', temp_file_name,
                '-d', taxdb,
            ],
            capture_output=True,
            text=True
        )

    taxonomy_data = {}
    for line in result.stdout.strip().split('\n'):
        fields = line.split('\t')[1:]
        if len(fields) == 3:
            taxid, taxon_details, ranks = fields[0], fields[1], fields[2]
            lineage_list = taxon_details.split(';')
            ranks_list = ranks.split(';')
            taxonomy = {
                rank: name for rank,
                name in zip(ranks_list, lineage_list)
                if rank in TAXONOMIC_RANKS
            }
            taxonomy_data[taxid] = taxonomy
        else:
            print("Warning: Unexpected format in taxonkit stdout."
                  f" This may result in missing taxonomy information:\n{line}",
                  file=sys.stderr)
    return taxonomy_data


if __name__ == '__main__':
    main()
