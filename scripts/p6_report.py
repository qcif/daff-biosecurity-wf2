"""Build the workflow report."""

import argparse

from src.report import report
from src.utils import existing_path
from src.utils.config import Config

config = Config()


def main():
    """Build the workflow report."""
    args = _parse_args()
    config.configure(args.output_dir)
    config.configure_query_logger(args.query_dir)
    report.render(args.query_dir)


def _parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "query_dir", type=existing_path, help="Path to query output directory")
    parser.add_argument(
        "--output_dir",
        type=existing_path,
        default=config.output_dir,
        help=f"Path to output directory. Defaults to {config.output_dir}.")
    return parser.parse_args()


if __name__ == '__main__':
    main()
