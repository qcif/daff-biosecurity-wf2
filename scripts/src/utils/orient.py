"""Translate a DNA sequence.

# Search for COX1 domain:
hmmsearch --tblout out.tbl --noali pf00115.hmm query.fasta

# Fields in out.tbl:

target name                forward_frame_1
accession                  -
query name                 COX1
accession                  PF00115.25

--- full sequence: ---
E-value                    4.8e-59
score                      188.9
bias                       21.7

--- best 1 domain: ---
E-value                    5.3e-59
score                      188.7
bias                       21.7

--- domain number estimation: ---
exp                        1.0
reg                        1
clu                        0
ov                         0
env                        1
dom                        1
rep                        1
inc                        1
description of target      -

"""

import logging
import os
import subprocess
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)

COX_PROFILE_PATH = Path(__file__).parent / "ref_data/pf00115.hmm"


def orientate(queries: list[SeqRecord]) -> list[SeqRecord]:
    """Orientate given sequences using HMMSearch.

    For sequences that can't be oriented with HMMSearch due to missing COX1
    domain, the function will return both strands (+/-). Each strand will
    have seq.annotations["oriented"]=False if orientation was not possible,
    and seq.annotations["reverse_complement"]=True to identify the reverse
    strand.
    """
    aa_frames = []
    oriented_seqs = []
    query_index = {
        seq.id: seq
        for seq in queries
    }
    for seq in query_index.values():
        frames = translate(seq.seq)
        for i, s in enumerate(frames):
            frame = i - 3
            if frame >= 0:
                frame += 1
            aa_frames.append(
                SeqRecord(
                    s,
                    id=seq.id,
                    description=seq.description,
                    annotations={
                        "frame": frame,
                        "forward": frame > 0,
                    },
                )
            )
    filtered_seqs = search_cox_profile(aa_frames)
    filtered_seqids = {s.id for s in filtered_seqs}
    unmatched_queries = [
        seq for seq in query_index.values()
        if seq.id not in filtered_seqids
    ]
    for aa_seq in filtered_seqs:
        qseq = query_index[aa_seq.id]
        forward = aa_seq.annotations["forward"]
        seq = qseq.seq if forward else qseq.seq.reverse_complement()
        oriented_seqs.append(
            SeqRecord(
                seq,
                id=qseq.id,
                description=qseq.description,
                annotations={
                    "frame": aa_seq.annotations["frame"],
                    "forward": forward,
                    "oriented": True,
                },
            )
        )
    for seq in unmatched_queries:
        seq.annotations["oriented"] = False
        seq.annotations["reverse_complement"] = False
        oriented_seqs.append(seq)
        reverse_comp = seq.reverse_complement(
            id=True,
            name=True,
            description=True,
            annotations=True,
        )
        reverse_comp.annotations["reverse_complement"] = True
        oriented_seqs.append(reverse_comp)

    return oriented_seqs


def translate(raw_seq: str) -> dict[str, list[str]]:
    """Translate a DNA sequence into six forward/reverse frames."""
    seq = Seq(raw_seq.upper().replace("-", "").replace("\n", ""))
    strands = [
        seq,
        seq.reverse_complement(),
    ]
    mito_tables = [
        CodonTable.ambiguous_dna_by_id[2],   # vertebrate mt
        CodonTable.ambiguous_dna_by_id[5],   # invertebrate mt
    ]
    translated = []
    for seq in strands:
        for frame in range(3):
            for table in mito_tables:
                aa = seq[frame:].translate(table, to_stop=False)
                if aa:
                    break
            translated.append(aa)

    return translated


def search_cox_profile(seqlist: list[SeqRecord]) -> list[SeqRecord]:
    """Search for COX1 domain using hmmsearch.

    The provided SeqRecords should have a frame=< -3:3 > annotation with all
    six translation frames for each sequence. The id should be identical for
    all six frames. The function will return a list of SeqRecords that match
    the COX1 domain profile. The SeqRecords will be filtered to only include
    those that match the COX1 domain profile.
    """
    filtered_seqs = []
    seq_index = {
        seq.id: seq
        for seq in seqlist
    }
    tmp_name = None
    tmp_out_name = None

    try:
        with (
            tempfile.NamedTemporaryFile(
                prefix="hmmsearch_",
                suffix=".fasta",
                mode='w',
                delete=False,
                delete_on_close=False,
            ) as tmp,
            tempfile.NamedTemporaryFile(
                mode='r',
                delete=False,
                prefix="hmmsearch_",
                suffix=".tbl",
            ) as tmp_out
        ):
            tmp_name = tmp.name
            SeqIO.write(seqlist, tmp, "fasta")
            tmp.close()
            tmp_out_name = tmp_out.name
            subprocess.run([
                    "hmmsearch",
                    "--tblout",
                    tmp_out.name,
                    "--noali",
                    COX_PROFILE_PATH,
                    tmp_name,
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            lines = tmp_out.readlines()
            for line in lines:
                if line.startswith("#"):
                    continue
                fields = line.split()
                if len(fields) > 0 and float(fields[4]) < 1e-5:
                    seq = seq_index.get(fields[0])
                    # TODO: check how hmmsearch handles query with no hit
                    if seq:
                        filtered_seqs.append(seq)
                    else:
                        raise RuntimeError(
                            f"HMMSearch: Sequence ID {fields[0]} was returned"
                            " by HMMSearch but not found in the input FASTA.")

    except subprocess.CalledProcessError as e:
        logger.error(
            "HMMSearch: Error running HMMSearch. Please check the input"
            " sequences and the HMM profile.")
        logger.error(e.stderr)
        raise

    finally:
        for path in (tmp_name, tmp_out_name):
            if path and os.path.exists(path):
                os.remove(path)

    missing_seqs = set(seq_index.keys()) - {seq.id for seq in filtered_seqs}
    if missing_seqs:
        logger.warning(
            "HMMSearch: The following sequences could not be oriented because"
            " they did not match the COX1 HMM profile. Both orientations will"
            " be submitted to BOLD to determine orientation:"
            f" {', '.join(missing_seqs)}")

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
