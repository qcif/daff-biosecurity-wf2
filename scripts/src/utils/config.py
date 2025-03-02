"""Runtime configuration for the workflow.

All configuration values can be overridden with environment variables.
"""

import os
import csv
import json
import logging
from Bio import SeqIO
from datetime import datetime
from logging.config import dictConfig
from pathlib import Path

from .log import get_logging_config
from .utils import path_safe_str
from . import countries

logger = logging.getLogger(__name__)

REPORT_FILENAME = "report_{sample_id}_{timestamp}.html"
QUERY_DIR_PREFIX = 'query_'


class Config:

    USER_EMAIL = os.getenv("USER_EMAIL")
    NCBI_API_KEY = os.getenv("NCBI_API_KEY")
    TIMESTAMP_FILENAME = os.getenv("TIMESTAMP_FILENAME", 'timestamp.txt')
    INPUT_FASTA_FILEPATH = Path(os.getenv("INPUT_FASTA_FILEPATH",
                                          "tests/test-data/query.fasta"))
    ACCESSIONS_FILENAME = os.getenv("ACCESSIONS_FILENAME", "accessions.txt")
    TAXONOMY_FILE = os.getenv("TAXONOMY_FILENAME", 'taxonomy.csv')
    QUERY_TITLE_FILE = os.getenv("QUERY_TITLE_FILENAME", 'query_title.txt')
    HITS_JSON = os.getenv("HITS_JSON_FILENAME", 'hits.json')
    HITS_FASTA = os.getenv("HITS_FASTA_FILENAME", 'hits.fasta')
    TAXONOMY_ID_CSV = os.getenv("TAXONOMY_ID_CSV_FILENAME",
                                'assigned_taxonomy.csv')
    CANDIDATES_FASTA = os.getenv("CANDIDATES_FASTA_FILENAME",
                                 'candidates.fasta')
    CANDIDATES_CSV = os.getenv("CANDIDATES_CSV_FILENAME", 'candidates.csv')
    CANDIDATES_JSON = os.getenv("CANDIDATES_JSON_FILENAME", 'candidates.json')
    CANDIDATES_COUNT_FILE = os.getenv("CANDIDATES_COUNT_FILENAME",
                                      'candidates_count.txt')
    CANDIDATES_SOURCES_JSON = os.getenv("CANDIDATES_SOURCES_JSON_FILENAME",
                                        'candidates_sources.json')
    TOI_DETECTED_CSV = os.getenv("TOI_DETECTED_CSV_FILENAME",
                                 'taxa_of_concern_detected.csv')
    PMI_MATCH_CSV = os.getenv("PMI_MATCH_CSV_FILENAME",
                              'preliminary_id_match.csv')
    BOXPLOT_IMG = os.getenv("BOXPLOT_IMG_FILENAME",
                            'boxplot_image.png')
    DB_COVERAGE_JSON = os.getenv("DB_COVERAGE_JSON_FILENAME",
                                 'db_coverage.json')
    DB_COVERAGE_TOI_LIMIT = int(os.getenv("DB_COVERAGE_TOI_LIMIT", 10))
    TAXONKIT_DATA = Path(
        os.getenv("TAXONKIT_DATA", '~/.taxonkit')
    ).expanduser().absolute()

    DB_COVERAGE_MAX_CANDIDATES = 3
    FLAG_FILE_TEMPLATE = 'flag_{identifier}.txt'
    GBIF_LIMIT_RECORDS = int(os.getenv("GBIF_LIMIT_RECORDS", 500))
    GBIF_ACCEPTED_STATUS = os.getenv(
        "GBIF_ACCEPTED_STATUS",
        'accepted,doubtful',
    ).upper().replace(' ', '').split(',')
    LOG_FILENAME = 'run.log'
    QUERY_LOG_FILENAME = 'query.log'
    ENTREZ_LOCK_FILE = 'entrez.lock'
    ENTREZ_MAX_RETRIES = 3
    GBIF_FAST_LOCK_FILE = 'gbif-fast.lock'
    GBIF_SLOW_LOCK_FILE = 'gbif-slow.lock'
    GBIF_MAX_RETRIES = 3
    ERRORS_DIR = 'errors'

    class INPUTS:
        METADATA_CSV_HEADER = {
            "sample_id": "sample_id",
            "locus": "locus",
            "preliminary_id": "preliminary_id",
            "taxa_of_interest": "taxa_of_interest",
            "country": "country",
            "host": "host",
        }
        FASTA_FILEPATH = Path(
            os.getenv(
                "INPUT_FASTA_PATH",
                Path(__file__).parent.parent.parent.parent
                / 'tests/test-data/queries.fasta')
        )
        METADATA_PATH = Path(
            os.getenv(
                "INPUT_METADATA_CSV_FILEPATH",
                Path(__file__).parent.parent.parent.parent
                / 'tests/test-data/metadata.csv')
        )

    class ALIGNMENT:
        MIN_NT = int(os.getenv('MIN_NT', 400))
        MIN_Q_COVERAGE = float(os.getenv('MIN_Q_COVERAGE', 0.85))
        MIN_IDENTITY = float(os.getenv('MIN_IDENTITY', 0.935))
        MIN_IDENTITY_STRICT = float(os.getenv('MIN_IDENTITY_STRICT', 0.985))

    class OUTPUTS:
        TOI_DETECTED_HEADER = [
            "Taxon of interest",
            "Match rank",
            "Match taxon",
            "Match species",
            "Match accession",
            "Match identity",
        ]

    class REPORT:
        TITLE = "Taxonomic assignment report"

    def configure(self, output_dir=None):
        if output_dir:
            self.set_output_dir(output_dir)
        conf = get_logging_config(self.output_dir / self.LOG_FILENAME)
        dictConfig(conf)

    @property
    def output_dir(self):
        return Path(os.getenv("OUTPUT_DIR", 'output'))

    def set_output_dir(self, output_dir):
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)
        os.environ["OUTPUT_DIR"] = str(output_dir)

    def configure_query_logger(self, query_dir):
        conf = get_logging_config(query_dir / self.QUERY_LOG_FILENAME)
        dictConfig(conf)

    def get_query_ix(self, ix_or_dir):
        """Resolve query index/dir to query index."""
        if (
            isinstance(ix_or_dir, str) and QUERY_DIR_PREFIX in ix_or_dir
        ) or isinstance(ix_or_dir, Path):
            query_dir = Path(ix_or_dir)
            return int(query_dir.name.split("_")[1]) - 1
        return ix_or_dir

    def get_query_dir(self, ix_or_dir):
        """Resolve query index/dir to query dir Path."""
        if (
            isinstance(ix_or_dir, str) and QUERY_DIR_PREFIX in ix_or_dir
        ) or isinstance(ix_or_dir, Path):
            query_dir = Path(ix_or_dir)
            return query_dir

        query_ix = int(ix_or_dir)
        sample_id = self.get_sample_id(query_ix)
        query_dir = (
            self.output_dir
            / f"{QUERY_DIR_PREFIX}{query_ix + 1:>03}_{sample_id}"
        )
        query_dir.mkdir(exist_ok=True, parents=True)
        return query_dir

    def get_sample_id(self, query):
        """Resolve query index/dir to sample ID."""
        query_ix = self.get_query_ix(query)
        return self.read_query_fasta(query_ix).id.split('.')[0]

    @property
    def taxonomy_path(self):
        return self.output_dir / self.TAXONOMY_FILE

    @property
    def entrez_lock_file(self):
        return self.output_dir / self.ENTREZ_LOCK_FILE

    @property
    def gbif_slow_lock_file(self):
        return self.output_dir / self.GBIF_SLOW_LOCK_FILE

    @property
    def gbif_fast_lock_file(self):
        return self.output_dir / self.GBIF_FAST_LOCK_FILE

    @property
    def start_time(self) -> datetime:
        path = self.output_dir / self.TIMESTAMP_FILENAME
        if path.exists():
            ts = path.read_text().strip(' \n')
            return datetime.strptime(ts, "%Y%m%d %H%M%S")
        now = datetime.now()
        timestamp = now.strftime("%Y%m%d %H%M%S")
        path.write_text(timestamp)
        return now

    @property
    def metadata(self) -> dict[str, dict]:
        """Read metadata from CSV file."""
        if getattr(self, '_metadata', None):
            return self._metadata
        self._metadata = {}
        with self.INPUTS.METADATA_PATH.open() as f:
            header = self.INPUTS.METADATA_CSV_HEADER
            for row in csv.DictReader(f):
                sample_id = row.pop(header["sample_id"]).split('.')[0]
                self._metadata[sample_id] = {
                    key: row[colname].strip()
                    for key, colname in header.items()
                    if key != "sample_id"
                }
        return self._metadata

    def _get_metadata_for_query(self, query, field) -> str:
        sample_id = self.get_sample_id(query)
        return self.metadata[sample_id][
            self.INPUTS.METADATA_CSV_HEADER[field]
        ]

    def get_locus_for_query(self, query) -> str:
        return self._get_metadata_for_query(query, "locus")

    def get_pmi_for_query(self, query) -> str:
        return self._get_metadata_for_query(query, "preliminary_id")

    def get_country_code_for_query(self, query) -> str:
        country = self._get_metadata_for_query(query, "country")
        return countries.get_code(country)

    def get_toi_list_for_query(self, query) -> list[str]:
        """Read taxa of interest from TOI file."""
        toi_field = self._get_metadata_for_query(query, "taxa_of_interest")
        return toi_field.split('|')

    def get_report_path(self, query):
        query_ix = self.get_query_ix(query)
        return self.get_query_dir(query_ix) / path_safe_str(
            REPORT_FILENAME.format(
                sample_id=self.get_sample_id(query_ix).replace('.', '_'),
                timestamp='NOW',  # self.timestamp, # ! TODO
            )
        )

    def read_query_fasta(self, index=None) -> list[SeqIO.SeqRecord]:
        """Read query FASTA file."""
        if not hasattr(self, "query_sequences"):
            self.query_sequences = list(
                SeqIO.parse(self.INPUT_FASTA_FILEPATH, "fasta"))
        if index is not None:
            return self.query_sequences[int(index)]
        return self.query_sequences

    def read_blast_hits_json(self, query):
        """Read BLAST hits from JSON file."""
        query_dir = self.get_query_dir(query)
        path = query_dir / self.HITS_JSON
        return self.read_json(path)

    def read_blast_hits_fasta(self, query) -> list[SeqIO.SeqRecord]:
        """Read BLAST hits from JSON file."""
        query_dir = self.get_query_dir(query)
        path = query_dir / self.HITS_FASTA
        return self.read_fasta(path)

    def read_taxonomy_file(self) -> dict[str, dict[str, str]]:
        """Read taxonomy from CSV file."""
        taxonomies = {}
        with self.taxonomy_path.open() as f:
            for row in csv.DictReader(f):
                taxonomies[row["accession"].split('.')[0]] = row
        return taxonomies

    def read_json(self, path):
        """Read JSON file."""
        with path.open() as f:
            return json.load(f)

    def read_fasta(self, path):
        """Read FASTA file."""
        return list(SeqIO.parse(path, "fasta"))

    def to_json(self) -> dict:
        """Serialize object to JSON-friendly dict."""
        return {
            key: value
            for key, value in self.__dict__.items()
            if key not in (
                'query_sequences',
            )
        }
