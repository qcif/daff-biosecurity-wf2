"""Runtime configuration for the workflow."""

from .env import getenv


class Config:

    def __init__(self, output_dir):
        self.INPUT_HITS_JSON = output_dir / getenv("HITS_JSON")
        self.INPUT_HITS_FASTA = output_dir / getenv("HITS_FASTA")
        self.OUTPUT_FLAGS_CSV = output_dir / getenv("FLAGS_CSV")
        self.OUTPUT_ID_FILE = output_dir / getenv("ID_FILE")
        self.OUTPUT_CANDIDATES_FASTA = output_dir / getenv("CANDIDATES_FASTA")
        self.OUTPUT_CANDIDATES_CSV = output_dir / getenv("CANDIDATES_CSV")
