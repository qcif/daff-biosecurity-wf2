"""Use the BOLD API to search for similar sequences to query."""

import argparse
import logging

from pathlib import Path

from src.bold.id_engine import BoldSearch
from src.utils import existing_path
from src.utils.config import Config

logger = logging.getLogger(__name__)
config = Config()


def main():
    args = _parse_args()
    logger.info(f"Searching BOLD with query {args.fasta_file}...")
    result = BoldSearch(args.fasta_file)
    print('Done')


def _parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "fasta_file",
        type=existing_path,
        help="Path to the FASTA file containing sequences to search.",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        help="Directory to save parsed output files (JSON and FASTA). Defaults"
             f" to env variable 'OUTPUT_DIR' or '{config.output_dir}'.",
        default=config.output_dir,
    )
    return parser.parse_args()


if __name__ == '__main__':
    main()
