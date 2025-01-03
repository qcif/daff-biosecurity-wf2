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
from Bio.Blast import NCBIXML


def parse_blast_xml(xml_file: str, output_file: str = None):
    """Parse BLAST XML output file and extract information about alignments.

    Args:
        xml_file (str): Path to the BLAST XML file.
        output_file (str): Optional path to save parsed output.
    """
    with open(xml_file, "r") as handle:
        blast_records = NCBIXML.parse(handle)

        results = []
        # Each record is for a different query - in practice this needs to work
        # with up to 200 query sequences
        for blast_record in blast_records:
            # each alignment is a "hit":
            for alignment in blast_record.alignments:

                # get alignment length by adding length of each hsp
                # get query coverage by dividing alignment length by query length

                # each hsp is a high-scoring pair (matching chunk of DNA):
                for hsp in alignment.hsps:
                    results.append({
                        "query": blast_record.query,
                        # use this to pull taxon info using blastdbcmd tool:
                        "subject": alignment.hit_def,
                        "score": hsp.score,
                        "e_value": hsp.expect,
                        "query_start": hsp.query_start,
                        "query_end": hsp.query_end,
                        "subject_start": hsp.sbjct_start,
                        "subject_end": hsp.sbjct_end,
                    })

        if output_file:
            with open(output_file, "w") as out_handle:
                for result in results:
                    out_handle.write(
                        f"Query: {result['query']}\n"
                        f"Subject: {result['subject']}\n"
                        f"Score: {result['score']}\n"
                        f"E-value: {result['e_value']}\n"
                        f"Query range: {result['query_start']}-{result['query_end']}\n"
                        f"Subject range: {result['subject_start']}-{result['subject_end']}\n\n"
                    )
        else:
            for result in results:
                print(
                    f"Query: {result['query']}\n"
                    f"Subject: {result['subject']}\n"
                    f"Score: {result['score']}\n"
                    f"E-value: {result['e_value']}\n"
                    f"Query range: {result['query_start']}-{result['query_end']}\n"
                    f"Subject range: {result['subject_start']}-{result['subject_end']}\n"
                )

def main():
    parser = argparse.ArgumentParser(description="Parse BLAST XML output file.")
    parser.add_argument(
        "xml_file",
        type=str,
        help="Path to the BLAST XML file to parse.",
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Optional path to save the parsed output to a file.",
    )

    args = parser.parse_args()

    parse_blast_xml(args.xml_file, args.output)

if __name__ == "__main__":
    main()
