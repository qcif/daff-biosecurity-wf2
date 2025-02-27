"""Analyze the database coverage of target species at the given locus.

Database coverage is analysed at three levels:

1. Target species coverage: The number of records for the target species
2. Related species coverage: The number of records for species related to the
   target species
3. Related species from sample country of origin: as for (2), but only for
   species which have occurence records in the same country as the target
   species.
"""

import argparse
import json
import logging

from src.utils import existing_path
from src.utils.config import Config
from src.coverage import assess_coverage

logger = logging.getLogger(__name__)
config = Config()

MODULE_NAME = "Database Coverage"


def main():
    args = parse_args()
    config.configure(args.output_dir)
    config.configure_query_logger(args.query_dir)
    results = assess_coverage(args.query_dir)
    write_db_coverage(args.query_dir, results)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "query_dir", type=existing_path, help="Path to query output directory")
    parser.add_argument(
        "--output_dir",
        type=existing_path,
        default=config.output_dir,
        help=f"Path to output directory. Defaults to {config.output_dir}.")
    return parser.parse_args()


def write_db_coverage(query_dir, results):
    path = query_dir / config.DB_COVERAGE_JSON
    with path.open("w") as f:
        json.dump(results, f, indent=2)
    logger.info(
        f"[{MODULE_NAME}]: Database coverage data written to {path}")
    return path


if __name__ == '__main__':
    main()
