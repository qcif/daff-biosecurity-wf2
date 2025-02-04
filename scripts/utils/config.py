"""Runtime configuration for the workflow."""

import os
import csv
import json
from Bio import SeqIO
from pathlib import Path
from .flags import FLAG_DETAILS


class Config:

    output_dir = Path(os.getenv("OUTPUT_DIR", 'output'))
    INPUT_TOI_FILEPATH = Path(
        os.getenv(
            "INPUT_TOI_FILEPATH",
            'tests/test-data/taxa_of_interest.txt')
    )
    ACCESSIONS_FILENAME = os.getenv("ACCESSIONS_FILENAME", "accessions.txt")
    TAXONOMY_FILE = os.getenv("TAXONOMY_FILENAME", 'taxonomy.csv')
    QUERY_TITLE_FILE = os.getenv("QUERY_TITLE_FILENAME", 'query_title.txt')
    HITS_JSON = os.getenv("HITS_JSON_FILENAME", 'hits.json')
    HITS_FASTA = os.getenv("HITS_FASTA_FILENAME", 'hits.fasta')
    FLAGS_CSV = os.getenv("FLAGS_CSV_FILENAME", 'flags.csv')
    TAXONOMY_ID_FILE = os.getenv("TAXONOMY_ID_FILENAME",
                                 'assigned_taxonomy.txt')
    CANDIDATES_FASTA = os.getenv("CANDIDATES_FASTA_FILENAME",
                                 'candidates.fasta')
    CANDIDATES_CSV = os.getenv("CANDIDATES_CSV_FILENAME", 'candidates.csv')
    CANDIDATES_JSON = os.getenv("CANDIDATES_JSON_FILENAME", 'candidates.json')
    TOI_DETECTED_CSV = os.getenv("TOI_DETECTED_CSV_FILENAME",
                                 'taxa_of_concern_detected.csv')

    class ALIGNMENT:
        MIN_NT = int(os.getenv('MIN_NT', 400))
        MIN_Q_COVERAGE = float(os.getenv('MIN_Q_COVERAGE', 0.85))
        MIN_IDENTITY = float(os.getenv('MIN_IDENTITY', 0.935))
        MIN_IDENTITY_STRICT = float(os.getenv('MIN_IDENTITY_STRICT', 0.985))

    class FLAGS:
        """Flags for reporting outcomes."""
        CANDIDATES = 1
        TOI = 2
        A = "A"
        B = "B"
        C = "C"
        D = "D"
        E = "E"

    @property
    def taxonomy_path(self):
        return self.output_dir / self.TAXONOMY_FILE

    def set_output_dir(self, output_dir):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True, parents=True)

    def read_blast_hits_json(self, query_dir):
        """Read BLAST hits from JSON file."""
        path = query_dir / self.HITS_JSON
        return self._read_json(path)

    def read_blast_hits_fasta(self, query_dir):
        """Read BLAST hits from JSON file."""
        path = query_dir / self.HITS_FASTA
        return self._read_fasta(path)

    def read_taxa_of_interest(self) -> list[str]:
        """Read taxa of interest from TOI file."""
        return [
            line.strip()
            for line in self.INPUT_TOI_FILEPATH.read_text().splitlines()
            if line.strip()
        ]

    def read_taxonomy_file(self):
        """Read taxonomy from CSV file."""
        taxonomies = {}
        with self.taxonomy_path.open() as f:
            for row in csv.DictReader(f):
                taxonomies[row["accession"].split('.')[0]] = row
        return taxonomies

    def write_flag(self, query_ix, flag_num, outcome):
        """Write flags to CSV file."""
        path = self.get_query_dir(query_ix) / self.FLAGS_CSV
        if not path.exists():
            with path.open("w") as f:
                writer = csv.writer(f)
                writer.writerow(["name", "number", "outcome", "explanation"])
        with path.open("a") as f:
            writer = csv.writer(f)
            writer.writerow([
                FLAG_DETAILS[flag_num]['name'],
                flag_num,
                outcome,
                FLAG_DETAILS[flag_num]['explanation'][outcome],
            ])

    def get_query_dir(self, query_ix):
        d = self.output_dir / f"query_{query_ix}"
        d.mkdir(exist_ok=True, parents=True)
        return d

    def ix_from_query_dir(self, query_dir):
        query_dir = Path(query_dir)
        return int(query_dir.name.split("_")[1])

    def _read_json(self, path):
        """Read JSON file."""
        with path.open() as f:
            return json.load(f)

    def _read_fasta(self, path):
        """Read FASTA file."""
        return list(SeqIO.parse(path, "fasta"))
