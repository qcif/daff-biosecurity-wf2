"""Run logic for pipeline phase 1-2.

- Attempt species ID from BLAST results.json
- Detect Taxon of concern?
- Write flags 1, 2
- Identify candidate species
- Write hits mapped to candidate species

"""

import argparse
import logging
from pathlib import Path
from types import SimpleNamespace
from utils.config import Config
from utils import deduplicate

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

MIN_IDENTITY_STRICT = 0.985
MIN_IDENTITY = 0.935
ALIGN_CRITERA = SimpleNamespace(
    nt=400,
    query_coverage=0.85,
)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "query_dir", type=Path, help="Path to query output directory")
args = parser.parse_args()
config = Config(args.output_dir)


def main():
    hits = config.read_blast_hits_json(args.query_dir)
    # sequences = config.read_blast_hits_fasta(args.query_dir)
    _assign_species_id(hits)
    _detect_taxa_of_interest(hits)


def _assign_species_id(hits):
    """Attempt species ID from BLAST hits.json data."""
    taxonomies = config.read_taxonomy_file()
    for query_ix, query in enumerate(hits):
        candidate_hits = [
            hit for hit in hits[query]
            if (
                hit["alignment_length"] >= ALIGN_CRITERA.nt
                or hit["query_coverage"] >= ALIGN_CRITERA.query_coverage
            ) and hit["identity"] >= MIN_IDENTITY
        ]
        for hit in candidate_hits:
            hit["species"] = _get_species(hit, taxonomies)
        candidate_hits_strict = [
            hit for hit in candidate_hits
            if hit["identity"] >= MIN_IDENTITY_STRICT
        ]
        candidate_species = deduplicate([
            hit for hit in candidate_hits
        ], key=lambda x: x["species"])
        candidate_species_strict = deduplicate([
            hit for hit in candidate_hits_strict
        ], key=lambda x: x["species"])
        _write_candidate_flags(
            query_ix,
            candidate_species_strict,
            candidate_species,
        )
        _write_candidates(
            query_ix,
            candidate_hits_strict or candidate_hits,
            candidate_species_strict or candidate_species,
        )


def _get_species(hit, taxonomies):
    """Retrieve species name for given hit."""
    tax = taxonomies.get(hit["accession"])
    if not tax:
        logger.warning(f"Taxonomy record not found for {hit['accession']} -"
                       " this hit could not be included in the candidate"
                       " species list.")
        return None
    return tax["species"]


def _detect_taxa_of_interest():
    """Cross-reference TOIs against candidates.

    Each TOI can be at any taxonomic level, so need to check at appropriate
    taxonomic level.
    """
    pass


def _write_candidate_flags(query, candidates_strict, candidates):
    if len(candidates_strict) == 1:
        flag = config.FLAGS.A
    elif len(candidates_strict) > 1 and len(candidates_strict) < 4:
        flag = config.FLAGS.B
    elif len(candidates_strict) >= 4:
        flag = config.FLAGS.C
    elif len(candidates):
        flag = config.FLAGS.D
    else:
        flag = config.FLAGS.E
    flag_number = 1
    config.write_flag(query, flag_number, flag)


def _write_candidates(
    query_ix: int,
    candidate_hits: list[dict],
    candidate_species: list[str],
):
    """Write candidates hits and species to file."""
    config.write_candidates(query_ix, candidate_hits, candidate_species)


def _write_toi_detected():
    """DOC."""
    pass


if __name__ == "__main__":
    main()
