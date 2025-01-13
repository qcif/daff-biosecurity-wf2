"""Parse BLAST XML output file.

This script parses a BLAST XML output file and writes the extracted data
to JSON and FASTA (DNA sequence) files.

JSON data should look something like this (* indicates value that may need to
be calculated):

records:
- query:
  length: int
  sequence: str (DNA sequence)
  title: str
  alignments:
  - *alignment_length: int
    *query_coverage: float
    *identity: float
    *score: float
    subject:
      *length: int
      *sequence: str (DNA sequence)
      title: str
      *taxon:
        *kingdom: str
        *phylum: str
        *class: str
        *order: str
        *family: str
        *genus: str
        *species: str
    hsps:
    - score: float
      e_value: float
      identity: float
      query_start: int
      query_end: int
      subject_start: int
      subject_end: int
      *length: int

"""

import argparse
import json
import os
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

TBD = "TBD"
# TODO: consider if the FORWARE_STRAND can be > 0 and REVERSE_STRAND can be < 0
FORWARD_STRAND = 1
REVERSE_STRAND = -1


def calculate_query_coverage(alignment_length, query_length):
    """Calculate query coverage as a percentage."""
    return (alignment_length / query_length) * 100 if query_length > 0 else 0


def parse_blast_xml(blast_xml_path: str, output_dir: str = None):
    """Parse BLAST XML output file and extract information about alignments.

    Args:
        xml_file (str): Path to the BLAST XML file.
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
            # each alignment is a "hit":
            for alignment in blast_record.alignments:
                query_coverage = calculate_query_coverage(
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
                    "query_coverage": query_coverage,
                    "score": TBD,
                    "e_value": TBD,
                    "identity_percent": TBD,
                    "hits": [],
                }

                # each hsp is a high-scoring pair (matching chunk of DNA):
                for hsp in alignment.hsps:
                    if hsp.query_start < hsp.query_end:
                        query_frame = FORWARD_STRAND
                    else:
                        query_frame = REVERSE_STRAND

                    if hsp.sbjct_start < hsp.sbjct_end:
                        hit_frame = FORWARD_STRAND
                    else:
                        hit_frame = REVERSE_STRAND
                    hsp_record = {
                        # use this to pull taxon info using blastdbcmd tool:
                        "score": hsp.score,
                        "e_value": hsp.expect,
                        "identity": hsp.identities,
                        "query_frame": query_frame,
                        "hit_frame": hit_frame,
                        "positive": hsp.positives,
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

        # Create output directory if it doesn't exist
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        blast_hits_json_path = os.path.join(output_dir, "blast_hits.json")
        with open(blast_hits_json_path, "w") as json_file:
            json.dump(results, json_file, indent=4)

        blast_hits_fasta_path = os.path.join(output_dir, "blast_hits.fasta")
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
        "output_dir",
        type=str,
        help="Directory to save parsed output files (JSON and FASTA).",
    )

    args = parser.parse_args()

    parse_blast_xml(args.blast_xml_path, args.output_dir)


if __name__ == "__main__":
    main()
