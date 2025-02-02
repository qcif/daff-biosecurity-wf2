"""Parse BLAST XML output file.
JSON data refer to blast_output_example.json.
"""

import argparse
import json
import logging
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pathlib import Path

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


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


def calculate_hit_score(hsps):
    """Calculate the total scores of all hsps for a hit."""
    return sum(hsp.score for hsp in hsps)


def calculate_hit_e_value(hit, effective_search_space):
    """Calculate the e_value for a hit."""
    if len(hit.hsps) == 1:
        return hit.hsps[0].expect
    return effective_search_space * 2 ** (-sum(hsp.bits for hsp in hit.hsps))


def calculate_hit_identity_percent(hsps):
    """Calculate the total identity of all hsps for a hit."""
    total_hsps_identity = sum(hsp.identities for hsp in hsps)
    total_hsp_align_length = sum((hsp.sbjct_end - hsp.sbjct_start)
                                 for hsp in hsps)
    hit_identity_percent = round(
        total_hsps_identity / total_hsp_align_length * 100,
        2,
    )
    return hit_identity_percent if total_hsp_align_length > 0 else 0


def calculate_hit_query_coverage_percent(alignment_length, query_length):
    """Calculate query coverage as a percentage."""
    return round(
        alignment_length / query_length * 100,
        2,
    ) if query_length > 0 else 0


def parse_blast_xml(
    blast_xml_path: str,
    output_dir: str = None,
):
    """Parse BLAST XML output file and extract information about alignments.

    Args:
        xml_file (str): Path to the BLAST XML file.
        output_file (str): Optional path to save parsed output.
    """
    with open(blast_xml_path, "r") as handle:
        blast_records = NCBIXML.parse(handle)
        results = []
        fasta_results = []

        for blast_record in blast_records:
            query_record = {
                "query_title": blast_record.query,
                "length": blast_record.query_length,
                "hits": []
            }
            effective_search_space = blast_record.effective_search_space

            for alignment in blast_record.alignments:
                hit_score = calculate_hit_score(alignment.hsps)
                hit_e_value = calculate_hit_e_value(
                    alignment,
                    effective_search_space)
                hit_identity_percent = calculate_hit_identity_percent(
                    alignment.hsps
                )
                hit_query_coverage = calculate_hit_query_coverage_percent(
                    alignment.length,
                    blast_record.query_length
                )
                hit_record = {
                    "hit_id": alignment.hit_id,
                    "hit_def": alignment.hit_def,
                    "accession": alignment.accession,
                    "subject_length": alignment.length,
                    "query_coverage_percent": hit_query_coverage,
                    "score": hit_score,
                    "e_value": hit_e_value,
                    "identity_percent": hit_identity_percent,
                    "hsps": [],
                }

                for hsp in alignment.hsps:
                    hsp_record = {
                        "score": hsp.score,
                        "e_value": hsp.expect,
                        "identity": hsp.identities,
                        "strand_query": hsp.strand[0],
                        "strand_subject": hsp.strand[1],
                        "gaps": hsp.gaps,
                        "query_start": hsp.query_start,
                        "query_end": hsp.query_end,
                        "subject_start": hsp.sbjct_start,
                        "subject_end": hsp.sbjct_end,
                        "alignment_length": hsp.align_length,
                        "alignment": {
                            "query": hsp.query,
                            "subject": hsp.sbjct,
                            "midline": hsp.match,
                        }
                    }
                    hit_record["hsps"].append(hsp_record)
                    fasta_results.append(SeqRecord(
                        Seq(hsp.sbjct),
                        id=alignment.accession,
                        description=alignment.hit_def))

                query_record["hits"].append(hit_record)
            results.append(query_record)

        all_accessions = list({
            hit["accession"]
            for query in results
            for hit in query["hits"]
        })

        output_dir.mkdir(exist_ok=True, parents=True)

        blast_hits_json_path = output_dir / "blast_hits.json"
        with open(blast_hits_json_path, "w") as f:
            json.dump(results, f, indent=4)
            logger.info(
                f"BLAST hits alignments written to {blast_hits_json_path}")

        blast_hits_fasta_path = output_dir / "blast_hits.fasta"
        with open(blast_hits_fasta_path, "w") as f:
            SeqIO.write(fasta_results, f, "fasta")
            logger.info(
                f"BLAST hits sequences written to {blast_hits_fasta_path}")

        hit_accesssions_path = output_dir / "accessions.txt"
        with open(hit_accesssions_path, "w") as f:
            f.write('\n'.join(all_accessions) + '\n')
            logger.info(
                f"BLAST hit accession IDs written to {hit_accesssions_path}")


def main():
    args = _parse_args()
    parse_blast_xml(args.blast_xml_path, args.output_dir)


if __name__ == "__main__":
    main()
