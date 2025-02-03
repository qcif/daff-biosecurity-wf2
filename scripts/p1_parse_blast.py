"""Parse BLAST output into per-query JSON and FASTA files."""

import argparse
import json
import logging
from Bio import SeqIO
from pathlib import Path

from blast.parse_xml import parse_blast_xml

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Number of accessions to collect per query for fetching taxonomy data:
COLLECT_N_ACCESSIONS = 10
BLAST_HITS_ACCESSIONS_FILENAME = "accessions.txt"
BLAST_HITS_FILENAME = "blast_hits.json"
BLAST_FASTA_FILENAME = "blast_hits.fasta"


def main():
    args = _parse_args()
    args.output_dir.mkdir(exist_ok=True, parents=True)
    hits, fastas = parse_blast_xml(args.blast_xml_path)
    _write_hits(hits, args.output_dir)
    _write_fastas(fastas, args.output_dir)
    _write_accessions(hits, args.output_dir)


def _parse_args():
    parser = argparse.ArgumentParser(
        description="Parse BLAST XML output file."
    )
    parser.add_argument(
        "blast_xml_path",
        type=str,
        help="Path to the BLAST XML file to parse.",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        help="Directory to save parsed output files (JSON and FASTA).",
        default=".",
    )
    return parser.parse_args()


def _get_query_dirname(i):
    return f"query_{i + 1}"


def _create_query_dir(query, output_dir, query_dirname):
    """Create a directory for this query and write the query title to file."""
    query_path = output_dir / query_dirname / 'query_title.txt'
    query_path.parent.mkdir(exist_ok=True, parents=True)
    with query_path.open("w") as f:
        f.write(query['query_title'])
        logger.info(f"BLAST query title written to {query_path}")


def _write_hits(hits, output_dir):
    """Write a JSON file of BLAST hits for each query sequence."""
    for i, query in enumerate(hits):
        query_dirname = _get_query_dirname(i)
        _create_query_dir(query, output_dir, query_dirname)
        path = output_dir / query_dirname / BLAST_HITS_FILENAME
        with path.open("w") as f:
            json.dump(hits, f, indent=2)
            logger.info(f"BLAST hits for query [{i}] written to {path}")


def _write_fastas(query_fastas, output_dir):
    """Write a fasta file of hit subjects for each query sequence."""
    for i, fastas in enumerate(query_fastas):
        if not fastas:
            continue
        query_dirname = _get_query_dirname(i)
        path = output_dir / query_dirname / BLAST_FASTA_FILENAME
        with open(path, "w") as f:
            SeqIO.write(fastas, f, "fasta")
            logger.info(
                f"BLAST hit sequences for query [{i}] written to {path}")


def _write_accessions(hits, output_dir):
    """Write a unique list of BLAST hit accession IDs to a file.

    These will be used for extracting taxonomy data.
    """
    hit_accesssions_path = output_dir / BLAST_HITS_ACCESSIONS_FILENAME
    all_accessions = list({
        hit["accession"]
        for query in hits
        for hit in query["hits"][:COLLECT_N_ACCESSIONS]
    })
    with open(hit_accesssions_path, "w") as f:
        f.write('\n'.join(all_accessions) + '\n')
        logger.info(
            f"BLAST hit accession IDs written to {hit_accesssions_path}")


if __name__ == "__main__":
    main()
