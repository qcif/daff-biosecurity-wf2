"""Runtime configuration for the workflow."""

import json
from Bio import SeqIO
from .env import getenv


class Config:

    def __init__(self, output_dir):
        self.INPUT_HITS_JSON = output_dir / getenv("HITS_JSON")
        self.INPUT_HITS_FASTA = output_dir / getenv("HITS_FASTA")
        self.OUTPUT_FLAGS_CSV = output_dir / getenv("FLAGS_CSV")
        self.OUTPUT_ID_FILE = output_dir / getenv("ID_FILE")
        self.OUTPUT_CANDIDATES_FASTA = output_dir / getenv("CANDIDATES_FASTA")
        self.OUTPUT_CANDIDATES_CSV = output_dir / getenv("CANDIDATES_CSV")

    def read_blast_hits_json(self):
        """Read BLAST hits from JSON file."""
        return self._read_json(self.INPUT_HITS_JSON)

    def read_blast_hits_fasta(self):
        """Read BLAST hits from JSON file."""
        return self._read_fasta(self.INPUT_HITS_FASTA)

    def _read_json(self, path):
        """Read JSON file."""
        with path.open() as f:
            return json.load(f)

    def _read_fasta(self, path):
        """Read FASTA file."""
        return list(SeqIO.parse(path, "fasta"))
