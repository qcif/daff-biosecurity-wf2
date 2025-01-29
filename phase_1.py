"""Run logic for pipeline phase 1-3.

- Attempt species ID from BLAST results.json
- Detect Taxon of concern?
- Write flags 1, 2
- Identify candidate species
- Write hits mapped to candidate species

"""

import argparse
from pathlib import Path
from utils.config import Config

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "toc_file", type=Path, help="Path to taxon of concern list file")
parser.add_argument(
    "output_dir", type=Path, help="Path to output directory")
args = parser.parse_args()
config = Config(args.output_dir)


def main():
    """Run phase 1-3 logic."""
    pass


def _assign_species_id():
    """Attempt species ID from BLAST results.json."""
    pass


def _detect_taxa_of_concern():
    """Attempt species ID from BLAST results.json."""
    pass


def _write_flags():
    """Attempt species ID from BLAST results.json."""
    pass


def _write_candidates():
    """Attempt species ID from BLAST results.json."""
    pass


if __name__ == "__main__":
    main()
