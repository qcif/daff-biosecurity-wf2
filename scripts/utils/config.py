"""Runtime configuration for the workflow."""

import csv
import json
from Bio import SeqIO
from .env import getenv

CANDIDATE_CSV_HEADER = [
    "species",
    "taxid",
    "accession",
    "hit_subject",  # rename
    "identity",  # rename, recalc
    "query_coverage",    # rename, recalc
    "alignment_length",  # TODO
    "e_value",
    "bitscore",  # rename
]


class Config:

    class FLAGS:
        """Flags for reporting outcomes."""
        A = "A"
        B = "B"
        C = "C"
        D = "D"
        E = "E"

    def __init__(self, output_dir):
        self.output_dir = output_dir
        self.INPUT_TOI_FILE = getenv("INPUT_TOI_FILE")
        self.HITS_JSON = getenv("HITS_JSON")
        self.HITS_FASTA = getenv("HITS_FASTA")
        self.TAXONOMY_FILE = output_dir / getenv("TAXONOMY")
        self.FLAGS_CSV = getenv("FLAGS_CSV")
        self.ID_FILE = getenv("ID_FILE")
        self.CANDIDATES_FASTA = getenv("CANDIDATES_FASTA")
        self.CANDIDATES_CSV = getenv("CANDIDATES_CSV")
        self.CANDIDATES_JSON = getenv("CANDIDATES_JSON")

    def read_blast_hits_json(self, query_dir):
        """Read BLAST hits from JSON file."""
        path = query_dir / self.HITS_JSON
        return self._read_json(path)

    def read_blast_hits_fasta(self, query_dir):
        """Read BLAST hits from JSON file."""
        path = query_dir / self.HITS_FASTA
        return self._read_fasta(path)

    def read_toi_file(self):
        """Read taxa of interest from TOI file."""
        return [
            line.strip()
            for line in self.INPUT_TOI_FILE.read_text().splitlines()
            if line.strip()
        ]

    def read_taxonomy_file(self):
        """Read taxonomy from CSV file."""
        taxonomies = {}
        with self.TAXONOMY_FILE.open() as f:
            for row in csv.DictReader(f):
                taxonomies[row["accession"]] = row
        return taxonomies

    def write_candidates(self, query_ix, hits, species):
        """Write candidate hits and species to CSV file."""
        self._write_candidates_json(query_ix, hits, species)
        self._write_candidates_csv(query_ix, hits, species)
        self._write_candidates_fasta(query_ix, species)

    def _write_candidates_json(self, query_ix, hits, species):
        path = (
            self.output_dir
            / self.get_query_dir(query_ix)
            / self.CANDIDATES_JSON)
        with path.open("w") as f:
            json.dump({
                "hits": hits,
                "species": species,
            }, f)

    def _write_candidates_csv(self, query_ix, species):
        path = (
            self.output_dir
            / self.get_query_dir(query_ix)
            / self.CANDIDATES_CSV)
        with path.open("w") as f:
            writer = csv.writer(f)
            writer.writerow(CANDIDATE_CSV_HEADER)
            for hit in species:
                writer.writerow([
                    hit.get(key, "")
                    for key in CANDIDATE_CSV_HEADER
                ])

    def _write_candidates_fasta(query_ix, species):
        """Write FASTA sequences for each candidate species to file."""
        pass

    def get_query_dir(self, query_ix):
        return self.output_dir / f"query_{query_ix}"

    def _read_json(self, path):
        """Read JSON file."""
        with path.open() as f:
            return json.load(f)

    def _read_fasta(self, path):
        """Read FASTA file."""
        return list(SeqIO.parse(path, "fasta"))
