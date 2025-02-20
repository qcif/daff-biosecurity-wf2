"""Analyze the diversity of reference sequence sources.

A source is defined as a publication or set of authors that are linked to the
genbank record for that sequence. If there are no references, no sources are
returned and the sequence is classified as "anonymous".

Many anonymous records are from automated genome annotation projects, often
carried out by NCBI themselves. These records are flagged so that the user can
be aware of the potential reduced credibility of these annotation.
"""

import argparse
import json
import logging

from src.entrez import genbank
from src.utils import existing_path, serialize
from src.utils.config import Config

logger = logging.getLogger(__name__)
config = Config()


def main():
    args = _parse_args()
    config.set_output_dir(args.output_dir)
    config.configure_query_logger(args.query_dir)
    species, hits = _read_candidate_hits(args.query_dir)
    species, hits = _collect_sources_per_species(species, hits)
    candidates = {
        "species": species,
        "hits": hits,
    }
    path = args.query_dir / config.CANDIDATES_JSON
    with path.open('w') as f:
        json.dump(candidates, f, default=serialize, indent=2)
    logger.info(f"Candidate hits with source diversity data written to {path}")


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


def _collect_sources_per_species(species, hits) -> list[dict]:
    accession_sources = genbank.fetch_sources([
        hit["accession"] for hit in hits
    ])
    for sp_ix, spec in enumerate(species):
        species_str = spec["species"]
        independent_sources = []
        for hit_ix, hit in enumerate(hits):
            if hit["species"] == species_str:
                source = accession_sources[hit["accession"]]
                if not any(
                    source.matches(agg_src)
                    for agg_src in independent_sources
                ):
                    independent_sources.append(source)
                hits[hit_ix]["source"] = source
        species[sp_ix]['independent_sources'] = len(independent_sources)
        species[sp_ix]['hit_count'] = hit_ix + 1
    return species, hits


def _read_candidate_hits(query_dir):
    candidates = config.read_json(query_dir / config.CANDIDATES_JSON)
    species = candidates["species"]
    hits = candidates["hits"]
    return species, hits


if __name__ == '__main__':
    main()
