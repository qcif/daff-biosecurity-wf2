"""Run logic for pipeline phase 1-2.

- Attempt species ID from BLAST results.json
- Detect Taxon of concern?
- Write flags 1, 2
- Identify candidate species
- Write hits mapped to candidate species

"""

import argparse
from pathlib import Path
from types import SimpleNamespace
from utils.config import Config

MIN_IDENTITY_STRICT = 0.985
MIN_IDENTITY = 0.95
ALIGN_CRITERA = SimpleNamespace(
    nt=400,
    query_coverage=0.85,
)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "toc_file", type=Path, help="Path to taxon of concern list file")
parser.add_argument(
    "output_dir", type=Path, help="Path to output directory")
args = parser.parse_args()
config = Config(args.output_dir)


def main():
    """Run phase 1-3 logic."""
    hits = config.read_blast_hits_json()
    sequences = config.read_blast_hits_fasta()
    _assign_species_id(hits)


def _assign_species_id(hits):
    """Attempt species ID from BLAST hits.json data."""
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
