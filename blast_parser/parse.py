"""Parse BLAST XML output file.
JSON data refer to blast_output_example.json.
"""

import argparse
import json
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pathlib import Path

from ncbi_taxonomy import NCBITaxonomy


def calculate_hit_score(hsps):
    """Calculate the total scores of all hsps for a hit."""
    return sum(hsp.score for hsp in hsps)


def calculate_hit_e_value(hit, effective_search_space):
    """Calculate the e_value for a hit."""
    if len(hit.hsps) == 1:
        return hit.hsps[0].expect
    return effective_search_space * 2 ** (-sum(hsp.bits for hsp in hit.hsps))


def calculate_hit_identity_percent(hsps, alignment_length):
    """Calculate the total identity of all hsps for a hit."""
    total_hsps_identity = sum(hsp.identities for hsp in hsps)
    total_hsp_align_length = sum(hsp.length for hsp in hsps)
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


def parse_blast_xml(blast_xml_path: str, input_db: str, output_dir: str = None):
    """Parse BLAST XML output file and extract information about alignments.

    Args:
        xml_file (str): Path to the BLAST XML file.
        input_db (str): Database path to use for retrieving taxon ID.
        output_file (str): Optional path to save parsed output.
    """
    with open(blast_xml_path, "r") as handle:
        blast_records = NCBIXML.parse(handle)

        results = []
        fasta_results = []
        # Each record is for a different query - in practice this needs to work
        # with up to 200 query sequences
        for blast_record in blast_records:
            query_record = {
                "query_title": blast_record.query,
                "length": blast_record.query_length,
                "hits": []
            }

            effective_search_space = blast_record.effective_search_space
            # each alignment is a "hit":
            for alignment in blast_record.alignments:
                hit_score = calculate_hit_score(alignment.hsps)
                hit_e_value = calculate_hit_e_value(
                    alignment,
                    effective_search_space)
                hit_identity_percent = calculate_hit_identity_percent(
                    alignment.hsps,
                    alignment.length,
                )
                hit_query_coverage = calculate_hit_query_coverage_percent(
                    alignment.length,
                    blast_record.query_length
                )
                # get alignment length by adding length of each hsp
                # get query coverage by dividing alignment length
                # by query length
                hit_record = {
                    "hit_id": alignment.hit_id,
                    "hit_def": alignment.hit_def,
                    "hit_accession": alignment.accession,
                    "subject_length": alignment.length,
                    "query_coverage_percent": hit_query_coverage,
                    "score": hit_score,
                    "e_value": hit_e_value,
                    "identity_percent": hit_identity_percent,
                    "hits": [],
                }

                # each hsp is a high-scoring pair (matching chunk of DNA):
                for hsp in alignment.hsps:
                    hsp_record = {
                        # use this to pull taxon info using blastdbcmd tool:
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
                    hit_record["hits"].append(hsp_record)

                    fasta_results.append(SeqRecord(
                        Seq(hsp.sbjct),
                        id=alignment.accession,
                        description=alignment.hit_def))
                query_record["hits"].append(hit_record)
            results.append(query_record)

        all_accessions = list({
            hit["hit_accession"]
            for query in results
            for hit in query["hits"]
        })
        taxonomies = NCBITaxonomy.extract(input_db, all_accessions)
        for query in results:
            for hit in query["hits"]:
                if hit["hit_accession"] in taxonomies:
                    hit["taxonomy"] = taxonomies[hit["hit_accession"]].json()

        output_dir.mkdir(exist_ok=True, parents=True)

        blast_hits_json_path = output_dir / "blast_hits.json"
        with open(blast_hits_json_path, "w") as json_file:
            json.dump(results, json_file, indent=4)

        blast_hits_fasta_path = output_dir / "blast_hits.fasta"
        with open(blast_hits_fasta_path, "w") as fasta_file:
            SeqIO.write(fasta_results, fasta_file, "fasta")


def main():
    parser = argparse.ArgumentParser(
        description="Parse BLAST XML output file."
    )
    parser.add_argument(
        "blast_xml_path",
        type=str,
        help="Path to the BLAST XML file to parse.",
    )
    parser.add_argument(
        "--input-db",
        type=Path,
        help="Database path to use for retrieving taxon ID.",
        default="input",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        help="Directory to save parsed output files (JSON and FASTA).",
        default="output",
    )

    args = parser.parse_args()
    parse_blast_xml(args.blast_xml_path, args.output_dir)


if __name__ == "__main__":
    main()
