"""Runtime configuration for the workflow.

All configuration values can be overridden with environment variables.
"""

import csv
import json
import logging
import os
import shutil
import tempfile
from datetime import datetime, timedelta
from logging.config import dictConfig
from pathlib import Path

from Bio import SeqIO

from . import countries
from .log import get_logging_config
from .utils import path_safe_str

logger = logging.getLogger(__name__)

MAP_FILENAME_TEMPLATE = "map_{taxon_str}.png"
REPORT_FILENAME = "report_{sample_id}_{timestamp}.html"
QUERY_DIR_PREFIX = 'query_'


class Config:

    USER_EMAIL = os.getenv("USER_EMAIL")
    NCBI_API_KEY = os.getenv("NCBI_API_KEY")
    TAXONKIT_DATA = os.getenv("TAXONKIT_DATA",
                              Path('~/.taxonkit').expanduser())
    TIMESTAMP_FILENAME = os.getenv("TIMESTAMP_FILENAME", 'timestamp.txt')
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
    INDEPENDENT_SOURCES_JSON = os.getenv("INDEPENDENT_SOURCES_JSON_FILENAME",
                                         'aggregated_sources.json')
    TOI_DETECTED_CSV = os.getenv("TOI_DETECTED_CSV_FILENAME",
                                 'taxa_of_concern_detected.csv')
    PMI_MATCH_CSV = os.getenv("PMI_MATCH_CSV_FILENAME",
                              'preliminary_id_match.csv')
    BOXPLOT_IMG_FILENAME = os.getenv("BOXPLOT_IMG_FILENAME",
                                     'identity-boxplot.png')
    TREE_NWK_FILENAME = os.getenv("TREE_NWK_FILENAME",
                                  'candidates.nwk')
    DB_COVERAGE_JSON = os.getenv("DB_COVERAGE_JSON_FILENAME",
                                 'db_coverage.json')
    DB_COVERAGE_TOI_LIMIT = int(os.getenv("DB_COVERAGE_TOI_LIMIT", 10))

    # BLAST-specific
    BLAST_MAX_TARGET_SEQS = int(os.getenv("BLAST_MAX_TARGET_SEQS", 2000))

    # BOLD-specific
    BOLD_TAXON_COUNT_JSON = os.getenv("BOLD_TAXON_COUNT_JSON",
                                      "bold_taxon_counts.json")
    BOLD_TAXON_COLLECTORS_JSON = os.getenv("BOLD_TAXON_COLLECTORS_JSON",
                                           "bold_taxon_collectors.json")
    BOLD_TAXONOMY_JSON = os.getenv("BOLD_TAXONOMY_JSON",
                                   "bold_taxonomy.json")

    BOLD_FLAG = 'BOLD'
    HMMSEARCH_MIN_EVALUE = 1e-5
    DB_COVERAGE_MAX_CANDIDATES = 3
    FLAG_FILE_TEMPLATE = '{identifier}.flag'
    GBIF_LIMIT_RECORDS = int(os.getenv("GBIF_LIMIT_RECORDS", 500))
    GBIF_ACCEPTED_STATUS = os.getenv(
        "GBIF_ACCEPTED_STATUS",
        'accepted,doubtful',
    ).upper().replace(' ', '').split(',')
    FLAG_DETAILS_CSV_PATH = (
        Path(__file__).parent.parent.parent.parent / 'flags.csv')
    LOG_FILENAME = 'run.log'
    QUERY_LOG_FILENAME = 'query.log'
    ENTREZ_CACHE_DIRNAME = 'entrez_cache'
    THROTTLE_SQLITE_FILE = 'throttle.sqlite'
    MAX_API_RETRIES = 3
    ERRORS_DIR = 'errors'
    TEMP_DIR_NAME = 'biosecurity'
    TEMP_CLEAN_AFTER_DAYS = 7

    TEMP_FILES = [
        ENTREZ_CACHE_DIRNAME,
    ]

    class INPUTS:
        METADATA_CSV_HEADER = {
            "sample_id": "sample_id",
            "locus": "locus",
            "preliminary_id": "preliminary_id",
            "taxa_of_interest": "taxa_of_interest",
            "country": "country",
            "host": "host",
        }
        METADATA_CSV_REQUIRED_FIELDS = (
            "sample_id",
            "preliminary_id",
        )
        FASTA_FILEPATH = Path(
            os.getenv(
                "INPUT_FASTA_FILEPATH",
                Path(__file__).parent.parent.parent.parent
                / 'tests/test-data/queries.fasta')
        )
        METADATA_PATH = Path(
            os.getenv(
                "INPUT_METADATA_CSV_FILEPATH",
                Path(__file__).parent.parent.parent.parent
                / 'tests/test-data/metadata.csv')
        )
        FASTA_MAX_LENGTH_NT = 3000
        FASTA_MIN_LENGTH_NT = 20
        FASTA_MAX_SEQUENCES = 150

    class CRITERIA:
        ALIGNMENT_MIN_NT = int(os.getenv('MIN_NT', 300))
        ALIGNMENT_MIN_Q_COVERAGE = float(os.getenv('MIN_Q_COVERAGE', 0.85))
        ALIGNMENT_MIN_IDENTITY = float(os.getenv('MIN_IDENTITY', 0.935))
        ALIGNMENT_MIN_IDENTITY_STRICT = float(
            os.getenv('MIN_IDENTITY_STRICT', 0.985))
        MAX_CANDIDATES_FOR_ANALYSIS = int(
            os.getenv('MAX_CANDIDATES_FOR_ANALYSIS', 3))
        SOURCES_MIN_COUNT = int(os.getenv('MIN_SOURCE_COUNT', 5))
        DB_COV_TARGET_MIN_A = int(os.getenv('DB_COV_MIN_A', 5))
        DB_COV_TARGET_MIN_B = int(os.getenv('DB_COV_MIN_B', 1))
        DB_COV_RELATED_MIN_A = int(os.getenv('DB_COV_RELATED_MIN_A', 90))
        DB_COV_RELATED_MIN_B = int(os.getenv('DB_COV_RELATED_MIN_B', 10))
        DB_COV_COUNTRY_MISSING_A = int(
            os.getenv('DB_COV_COUNTRY_MISSING_A', 1))

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
        TITLE = "Taxonomic identification report"

    def configure(self, output_dir=None, bold=False):
        if output_dir:
            self.set_output_dir(output_dir)
        conf = get_logging_config(self.output_dir / self.LOG_FILENAME)
        dictConfig(conf)
        if bold:
            self.bold_flag_file.write_text('1')

    @property
    def bold_flag_file(self):
        """Path to the BOLD flag file."""
        return self.output_dir / self.BOLD_FLAG

    @property
    def is_bold(self):
        """Check if this is a BOLD run."""
        return self.bold_flag_file.exists()

    @property
    def output_dir(self):
        return Path(os.getenv("OUTPUT_DIR", 'output'))

    def set_output_dir(self, output_dir):
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)
        os.environ["OUTPUT_DIR"] = str(output_dir)

    def create_query_dir(self, query_ix, query_title):
        """Create a directory for this query and write query title file."""
        query_dir = self.get_query_dir(query_ix)
        query_title_path = query_dir / self.QUERY_TITLE_FILE
        with query_title_path.open("w") as f:
            f.write(query_title)
            logger.info(f"Query title written to {query_title_path}")
        return query_dir

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
    def allowed_loci(self):
        path = (
            Path(os.getenv('ALLOWED_LOCI_PATH'))
            if os.getenv('ALLOWED_LOCI_PATH') else None)
        if not path:
            return []
        return [
            x.strip()
            for x in Path(path).read_text().split('\n')
            if x.strip()
        ]

    @property
    def taxonomy_path(self):
        return self.output_dir / self.TAXONOMY_FILE

    @property
    def entrez_cache_dir(self):
        return self.output_dir / self.ENTREZ_CACHE_DIRNAME

    @property
    def tempdir(self):
        tempdir = Path(tempfile.gettempdir()) / self.TEMP_DIR_NAME
        tempdir.mkdir(exist_ok=True, parents=True)
        return tempdir

    @property
    def user_tempdir(self):
        user_dir = self.tempdir / (self.USER_EMAIL or 'ANONYMOUS')
        user_dir.mkdir(exist_ok=True, parents=True)
        return user_dir

    @property
    def throttle_sqlite_path(self):
        return self.user_tempdir / self.THROTTLE_SQLITE_FILE

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
                    key: (
                        [
                            x.strip()
                            for x in row[colname].strip().split('|')
                            if x.strip()
                        ]
                        if 'interest' in key.lower()
                        else row[colname].strip()
                    )
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

    def get_country_for_query(self, query, code=False) -> str:
        country = self._get_metadata_for_query(query, "country")
        if not country:
            return None
        if code:
            return countries.get_code(country)
        return country

    def get_toi_list_for_query(self, query) -> list[str]:
        """Read taxa of interest from TOI file."""
        return self._get_metadata_for_query(query, "taxa_of_interest")

    def get_report_path(self, query):
        query_ix = self.get_query_ix(query)
        return self.get_query_dir(query_ix) / path_safe_str(
            REPORT_FILENAME.format(
                sample_id=self.get_sample_id(query_ix).replace('.', '_'),
                timestamp='NOW',  # self.timestamp, # ! TODO
            )
        )

    def read_flag_details_csv(self) -> dict[str, dict[str, dict]]:
        data = {}
        with self.FLAG_DETAILS_CSV_PATH.open() as f:
            reader = csv.DictReader(f)
            for row in reader:
                flag = data.get(row['id'], {'name': row['name']})
                flag['explanation'] = flag.get('explanation', {})
                flag['outcome'] = flag.get('outcome', {})
                flag['level'] = flag.get('level', {})
                flag['explanation'][row['value']] = row['explanation']
                flag['outcome'][row['value']] = row['outcome']
                flag['level'][row['value']] = int(row['level'])
                data[row['id']] = flag
        return data

    def read_query_fasta(self, index=None) -> list[SeqIO.SeqRecord]:
        """Read query FASTA file."""
        if not hasattr(self, "query_sequences"):
            self.query_sequences = list(
                SeqIO.parse(self.INPUTS.FASTA_FILEPATH, "fasta"))
        if index is not None:
            return self.query_sequences[int(index)]
        return self.query_sequences

    def read_hits_json(self, query):
        """Read BLAST hits from JSON file."""
        query_dir = self.get_query_dir(query)
        path = query_dir / self.HITS_JSON
        return self.read_json(path)

    def read_hits_fasta(self, query) -> list[SeqIO.SeqRecord]:
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

    def url_from_accession(self, accession):
        return f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}"

    def get_map_filename_for_target(self, taxon):
        return MAP_FILENAME_TEMPLATE.format(taxon_str=path_safe_str(taxon))

    def cleanup(self):
        """Remove temporary files."""
        logger.info("Cleaning temporary files from output dir...")
        for filename in self.TEMP_FILES:
            path = self.output_dir / filename
            if path.exists():
                if path.is_dir():
                    shutil.rmtree(path)
                else:
                    path.unlink()
        for path in self.tempdir.glob('*'):
            if path.is_dir():
                mtime = get_latest_mtime(path)
                if (
                    mtime
                    < datetime.now()
                    - timedelta(days=self.TEMP_CLEAN_AFTER_DAYS)
                ):
                    shutil.rmtree(path)


def get_latest_mtime(path: str) -> datetime:
    """Return the latest modification time of files/dirs within given dir."""
    latest_mtime = os.path.getmtime(path)
    for root, dirs, files in os.walk(path):
        for name in dirs + files:
            full_path = os.path.join(root, name)
            try:
                mtime = os.path.getmtime(full_path)
                if mtime > latest_mtime:
                    latest_mtime = mtime
            except FileNotFoundError:
                # Skip files that may have been deleted during walk
                continue

    return datetime.fromtimestamp(latest_mtime)
