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
import sys

from src.utils import existing_path, errors
from src.utils.config import Config
from src.coverage import assess_coverage

logger = logging.getLogger(__name__)
config = Config()

MODULE_NAME = "Database Coverage"


def main():
    args = parse_args()
    config.configure(args.output_dir, query_dir=args.query_dir)
    results, error_detected = assess_coverage(
        args.query_dir,
        is_bold=args.bold,
    )
    write_db_coverage(args.query_dir, results)
    config.cleanup()
    if error_detected:
        error_file_path = args.query_dir / 'error.p5.log'
        sys.stderr.write(
            f'[Query {args.query_dir.name}] An error occurred during database'
            ' coverage assessment that'
            ' prevented one or more target species from being assessed.'
            ' For further details, please consult the workflow report or error'
            f' file: {error_file_path}\n'
        )
        errors.report(
            error_file_path,
            query_dir=args.query_dir,
            location_min=errors.LOCATIONS.DATABASE_COVERAGE,
            location_max=errors.LOCATIONS.DB_COVERAGE_RELATED_COUNTRY,
        )


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "query_dir", type=existing_path, help="Path to query output directory")
    parser.add_argument(
        "--output_dir",
        type=existing_path,
        default=config.output_dir,
        help=f"Path to output directory. Defaults to {config.output_dir}.")
    parser.add_argument(
        "--bold",
        action="store_true",
        help="Reference the BOLD database instead of GenBank.")
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
