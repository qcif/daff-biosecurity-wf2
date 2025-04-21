"""Translate a DNA sequence.

# Search for COX1 domain:
hmmsearch --tblout out.tbl --noali pf00115.hmm query.fasta

# Example out.tbl:

#                                                               --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name        accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
forward_frame_1      -          COX1                 PF00115.25   4.8e-59  188.9  21.7   5.3e-59  188.7  21.7   1.0   1   0   0   1   1   1   1 -
#
# Program:         hmmsearch
# Version:         3.4 (Aug 2023)
# Pipeline mode:   SEARCH
# Query file:      pf00115.hmm
# Target file:     query.fasta
# Option settings: hmmsearch --tblout out.tbl --noali pf00115.hmm query.fasta
# Current dir:     /hmm
# Date:            Mon Apr 21 00:50:00 2025
# [ok]

"""

import logging
import os
import subprocess
import tempfile

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)


def translate(raw_seq: str) -> dict[str, list[str]]:
    """Translate a DNA sequence into six forward/reverse frames."""
    seq = Seq(raw_seq.upper().replace("-", "").replace("\n", ""))
    strands = {
        'f': seq,
        'r': seq.reverse_complement()
    }
    mito_tables = [
        CodonTable.ambiguous_dna_by_id[2],   # vertebrate mt
        CodonTable.ambiguous_dna_by_id[5],   # invertebrate mt
    ]
    translated = {}
    for strand, seq in strands.items():
        for frame in range(3):
            for table in mito_tables:
                aa = seq[frame:].translate(table, to_stop=False)
                if aa:
                    translated.setdefault(strand, []).append(aa)
                    break
                translated.get(strand, []).append(aa)
    return translated


def search_cox_profile(seqlist: list[SeqRecord]) -> list[SeqRecord]:
    """Search for COX1 domain using hmmsearch."""
    seq_index = {
        seq.id: seq
        for seq in seqlist
    }

    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as tmp:
        tmp_name = tmp.name
        SeqIO.write(seqlist, tmp, "fasta")
        tmp.close()

    cox_profile_path = os.getenv('COX1_HMM_PROFILE')
    if not cox_profile_path:
        raise EnvironmentError("COX1_HMM_PROFILE must be set to run HMMSearch")

    with tempfile.NamedTemporaryFile(delete=False, suffix=".tbl") as tmp_out:
        tmp_out_name = tmp_out.name
        subprocess.run([
                "hmmsearch",
                "--tblout",
                tmp_out.name,
                "--noali",
                cox_profile_path,
                tmp_name,
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        lines = tmp_out.readlines()
        filtered_seqs = []
        for line in lines:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) > 0 and float(fields[4]) < 1e-5:
                seq = seq_index.get(fields[0])
                if seq:
                    filtered_seqs.append(seq)
                else:
                    logger.error(
                        f"HMMSearch: Sequence ID {fields[0]} was returned by"
                        " HMMSearch but not found in the input file.")

    missing_seqs = set(seq_index.keys()) - {seq.id for seq in filtered_seqs}
    if missing_seqs:
        logger.error(
            "HMMSearch: The following sequences could not be oriented because"
            " they did not match the COX1 HMM profile:"
            f" {', '.join(missing_seqs)}")

    os.remove(tmp_name)
    os.remove(tmp_out_name)

    return filtered_seqs


if __name__ == '__main__':
    query = """
GGTAAAAAAAATGAGTTTTTGGCTTTTGCCTCCTTCTTTTCTTCTTTTATTGGCTTCTGCTGGTGTTGAAAGGGGTGTTGGTACTGGGTGAACTATTTATCCTCCTTTGTCAAGTGGTATTGCTCATTCTGGTGGTTCTGTTGATCTTGCTATTTTTTCTTTACATATAGCTGGTGCTTCTTCTATTATAGCTTCGATTAACTTTATAACTACAATAATAAAAATGCGTGCTCCTGGTATTTCTTTTGATCGTCTTTCTCTTTTTGTTTGATCAATTTTTATTACTACTTTTCTTCTTTTGTTATCTTTACCTGTTTTGGCTGGGGCTATAACAATGCTTCTTACTGATCGTAATGTAAATACTACTTTTTTTGATCCTGCTGGTGGTGGTGATCCTATTTTGTTTCAGCATTTATTTTGGTTTTTTGGTCATCCTGAAGTTTATATTTTAATATTACCTGGTTTTGGTATGATTTCTCATGTTGTTTCTCATTATTCTGGTAAGAGAGAGCCTTTTGGTTATTTGGGTATGGTTTATGCGATGGTTGCTATAGGTATACTTGGTTTTCTTGTTTGAGCACATCATATGTTTACTGTAGGTA
"""
    res = translate(query)
    for strand, frames in res.items():
        print(f"{strand}:")
        for i, frame in enumerate(frames):
            print(f"Frame {i}: {frame}")
        print()
