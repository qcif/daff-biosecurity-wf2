"""Run logic for pipeline phase 1-2.

- Attempt species ID from BLAST results.json (flag 1)
- Detect Taxa of Interest (flag 2)

Taxa of Interest output has the following CSV fields:

Taxon of interest: The provided TOI that matched a candidate species
Match rank: The taxonomic rank of the match (TOI rank may be above species)
Match taxon: The taxon that matched the TOI
Match species: The species of the candidate match
Match accession: The NCBI accession of the candidate match

"""

import argparse
import csv
import json
import logging
from Bio import SeqIO

from utils import deduplicate, existing_path
from utils.config import Config
from utils.flags import Flag, FLAGS

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

CANDIDATE_CSV_HEADER = [
    "species",
    "taxid",
    "accession",
    "hit_subject",
    "identity",
    "query_coverage",
    "alignment_length",
    "e_value",
    "bitscore",
]

config = Config()


def main():
    args = _parse_args()
    config.set_output_dir(args.output_dir)
    query = config.read_blast_hits_json(args.query_dir)
    candidate_species = _assign_species_id(query['hits'], args.query_dir)
    _detect_taxa_of_interest(candidate_species, args.query_dir)


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


def _assign_species_id(hits, query_dir):
    """Attempt species ID from BLAST hits.json data."""
    query_ix = config.ix_from_query_dir(query_dir)
    taxonomies = config.read_taxonomy_file()
    candidate_hits = [
        hit for hit in hits
        if (
            hit["alignment_length"] >= config.ALIGNMENT.MIN_NT
            or hit["query_coverage"] >= config.ALIGNMENT.MIN_Q_COVERAGE
        ) and hit["identity"] >= config.ALIGNMENT.MIN_IDENTITY
    ]
    for hit in candidate_hits:
        tax = taxonomies.get(hit["accession"])
        if tax:
            hit['taxonomy'] = tax
            hit["species"] = tax.get('species')
            hit['taxid'] = tax.get('taxid')
        else:
            hit['taxonomy'] = None
            hit["species"] = None
            hit['taxid'] = None
            logger.warning(
                f"Taxonomy record not found for {hit['accession']} -"
                " this hit could not be included in the candidate"
                " species list.")
    candidate_hits_strict = [
        hit for hit in candidate_hits
        if hit["identity"] >= config.ALIGNMENT.MIN_IDENTITY_STRICT
    ]
    candidate_species = deduplicate([
        hit for hit in candidate_hits
    ], key=lambda x: x.get("species"))
    candidate_species_strict = deduplicate([
        hit for hit in candidate_hits_strict
    ], key=lambda x: x.get("species"))

    # Add hit count for each candidate species
    for species in candidate_species:
        species["hit_count"] = sum(
            1 for hit in candidate_hits
            if hit["species"] == species["species"]
        )
    for species in candidate_species_strict:
        species["hit_count"] = sum(
            1 for hit in candidate_hits_strict
            if hit["species"] == species["species"]
        )

    _write_taxonomic_id(query_ix, candidate_hits_strict)
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
    return candidate_species_strict or candidate_species


def _write_taxonomic_id(query_ix, candidate_hits_strict):
    if len(candidate_hits_strict) != 1:
        logger.info(f"Query {query_ix} - not writing {config.TAXONOMY_ID_CSV}:"
                    " no taxonomic identification could be made.")
    else:
        query_dir = config.get_query_dir(query_ix)
        src = query_dir / config.CANDIDATES_CSV
        dest = query_dir / config.TAXONOMY_ID_CSV
        dest.write_text(src.read_text())


def _write_candidate_flags(query_ix, candidates_strict, candidates):
    if len(candidates_strict) == 1:
        value = FLAGS.A
    elif len(candidates_strict) > 1 and len(candidates_strict) < 4:
        value = FLAGS.B
    elif len(candidates_strict) >= 4:
        value = FLAGS.C
    elif len(candidates):
        value = FLAGS.D
    else:
        value = FLAGS.E
    Flag.write(query_ix, FLAGS.POSITIVE_ID, value)


def _write_candidates(
    query_ix: int,
    candidate_hits: list[dict],
    candidate_species: list[str],
):
    """Write candidates hits and species to file."""
    _write_candidates_json(query_ix, candidate_hits, candidate_species)
    _write_candidates_csv(query_ix, candidate_species)
    _write_candidates_fasta(query_ix, candidate_species)


def _write_candidates_json(query_ix, hits, species):
    path = config.get_query_dir(query_ix) / config.CANDIDATES_JSON
    with path.open("w") as f:
        json.dump({
            "hits": hits,
            "species": species,
        }, f, indent=2)
    logger.info(f"Written candidate hits/species to {path}")


def _write_candidates_csv(query_ix, species):
    path = config.get_query_dir(query_ix) / config.CANDIDATES_CSV
    with path.open("w") as f:
        writer = csv.writer(f)
        writer.writerow(CANDIDATE_CSV_HEADER)
        for hit in species:
            writer.writerow([
                hit.get(key, "")
                for key in CANDIDATE_CSV_HEADER
            ])
    logger.info(f"Written candidate species to {path}")


def _write_candidates_fasta(query_ix, species):
    """Write FASTA sequences for each candidate species to file."""
    query_dir = config.get_query_dir(query_ix)
    path = query_dir / config.CANDIDATES_FASTA
    fastas = config.read_blast_hits_fasta(query_dir)
    species_accessions = [hit["accession"] for hit in species]
    candidate_fastas = [
        fasta for fasta in fastas
        if fasta.id in species_accessions
    ]
    with open(path, "w") as f:
        SeqIO.write(candidate_fastas, f, "fasta")
    logger.info(f"Written candidate FASTA to {path}")


def _detect_taxa_of_interest(candidate_species, query_dir):
    """Cross-reference Taxa of Interest (TOIs) against candidate species.

    Each TOI can be at any taxonomic level, so need to check at appropriate
    taxonomic level.
    """
    taxa_of_interest = config.read_taxa_of_interest()
    if not taxa_of_interest:
        logger.info("No taxa of interest provided - no output written.")
        return
    candidate_taxa_flattened = [
        {
            'rank': rank,
            'taxon': taxon,
            'species': hit['species'],
            'accession': hit['accession'],
            'identity': hit['identity'],
        }
        for hit in candidate_species
        for rank, taxon in hit["taxonomy"].items()
    ]
    write_flag = True

    with (query_dir / config.TOI_DETECTED_CSV).open("w") as f:
        writer = csv.writer(f)
        writer.writerow(config.OUTPUTS.TOI_DETECTED_HEADER)

    for toi in taxa_of_interest:
        try:
            detected_taxon = [
                tax for tax in candidate_taxa_flattened
                if toi.lower() == tax['taxon'].lower()
            ][0]
        except IndexError:
            detected_taxon = {}
        _write_toi_detected(
            query_dir,
            toi,
            detected_taxon,
            write_flag=write_flag,
        )
        if detected_taxon:
            write_flag = False
    logger.info("Writing taxa of interest detection to"
                f" {query_dir / config.TOI_DETECTED_CSV}")


def _write_toi_detected(query_dir, toi, detected, write_flag=True):
    """Record whether a given TOI was detected and set flag."""
    if write_flag:
        if detected:
            value = FLAGS.B
        else:
            value = FLAGS.A
        Flag.write(
            config.ix_from_query_dir(query_dir),
            FLAGS.TOI,
            value,
        )
    with (query_dir / config.TOI_DETECTED_CSV).open("a") as f:
        writer = csv.writer(f)
        writer.writerow([
            toi,
            detected.get('rank', ''),
            detected.get('taxon', ''),
            detected.get('species', ''),
            detected.get('accession', ''),
            detected.get('identity', ''),
        ])


if __name__ == "__main__":
    main()
