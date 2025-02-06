"""Runtime configuration for the workflow."""

import json
from Bio import SeqIO
from .env import getenv


class Config:

    def __init__(self, output_dir):
        self.output_dir = output_dir
        self.HITS_JSON = getenv("HITS_JSON")
        self.HITS_FASTA = getenv("HITS_FASTA")
        self.FLAGS_CSV = getenv("FLAGS_CSV")
        self.ID_FILE = getenv("ID_FILE")
        self.CANDIDATES_FASTA = getenv("CANDIDATES_FASTA")
        self.CANDIDATES_CSV = getenv("CANDIDATES_CSV")
        self.ENTREZ_EMAIL = getenv("ENTRZ_EMAIL")

    def read_blast_hits_json(self):
        """Read BLAST hits from JSON file."""
        path = self._path_for_query(self.HITS_JSON)
        return self._read_json(path)

    def read_blast_hits_fasta(self):
        """Read BLAST hits from JSON file."""
        path = self._path_for_query(self.HITS_FASTA)
        return self._read_fasta(path)

    def _path_for_query(self, path, query_ix):
        """Return path for given with query index."""
        return self.output_dir / f"query_{query_ix}" / path

    def _read_json(self, path):
        """Read JSON file."""
        with path.open() as f:
            return json.load(f)

    def _read_fasta(self, path):
        """Read FASTA file."""
        return list(SeqIO.parse(path, "fasta"))
