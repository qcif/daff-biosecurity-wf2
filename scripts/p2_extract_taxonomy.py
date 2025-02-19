"""Extract taxids and taxonomic information from NCBI databases.

This requires access to the NCBI taxdump files (configurable by CLI param).

"""

import argparse
import csv
import logging
from pathlib import Path

try:
    from utils import existing_path
    from utils.config import Config
    from taxonomy import extract
    from taxonomy.extract import TAXONKIT_DATA, TAXONOMIC_RANKS
except ImportError:
    # For unit tests only
    from .utils import existing_path
    from .utils.config import Config
    from .taxonomy import extract
    from .taxonomy.extract import TAXONKIT_DATA, TAXONOMIC_RANKS

logger = logging.getLogger(__name__)
config = Config()


def main():
    args = _parse_args()
    with args.taxids_csv.open() as taxids_file:
        accession_taxids = {
            row[0]: row[1]
            for row in csv.reader(taxids_file)
        }
    taxids = sorted(set(accession_taxids.values()))
    taxonomies = extract.taxonomies(taxids, taxdb=args.taxdb_path)
    _write_csv(taxonomies, accession_taxids, args.output_csv)


def _parse_args():
    default_output_csv = config.output_dir / config.TAXONOMY_FILE
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'taxids_csv',
        type=existing_path,
        help='CSV file with columns (accession,taxid) to extract taxonomy'
             ' information for.',
    )
    parser.add_argument(
        '--taxdb',
        dest='taxdb_path',
        type=existing_path,
        help='Path to directory containing NCBI taxdump files for taxonkit.'
             f' Defaults to {TAXONKIT_DATA}',
        default=Path(TAXONKIT_DATA),
    )
    parser.add_argument(
        '--output',
        dest='output_csv',
        type=Path,
        help='CSV file where taxonomy data will be written. Defaults to'
             f' {default_output_csv}',
        default=default_output_csv,
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
    logger.info(f"Taxonomy records written to {output_csv}")


if __name__ == '__main__':
    main()
