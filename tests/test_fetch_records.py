import json
import unittest
from pathlib import Path

from scripts.entrez.genbank import fetch_sources, fetch_gb_records

EXPECTED_SINGLE_SOURCE = None

DATA_DIR = Path(__file__).parent / 'test-data'
EXPECTED_RECORD_COUNT = 1686
EXPECT_RECORDS_JSON = DATA_DIR / 'genbank_expect_ids.json'
EXPECT_SINGLE_SOURCE_JSON = DATA_DIR / 'genbank_single_source.json'
EXPECT_MULTIPLE_SOURCES_JSON = DATA_DIR / 'genbank_multiple_sources.json'

ACCESSION_1 = "NM_001126"
ACCESSION_2 = "HQ621368"
DATABASE = "nuccore"
SINGLE_ACCESSION = [ACCESSION_1]
MULTIPLE_ACCESSIONS = [ACCESSION_1, ACCESSION_2]
LOCUS = 'COI'
TAXID = "9606"


class TestFetchRecords(unittest.TestCase):

    # @patch('fetch_gb_records.fetch_metadata')
    def test_fetch_single_source(self):
        result = fetch_sources(SINGLE_ACCESSION)
        with EXPECT_SINGLE_SOURCE_JSON.open() as f:
            expected = json.load(f)
        self.assertEqual(result, expected)

    def test_fetch_multiple_source(self):
        result = fetch_sources(MULTIPLE_ACCESSIONS)
        with EXPECT_MULTIPLE_SOURCES_JSON.open() as f:
            expected = json.load(f)
        self.assertEqual(result, expected)

    # @patch('fetch_gb_records.fetch_')
    def test_fetch_gb_records_count(self):
        result = fetch_gb_records(LOCUS, TAXID, True)
        self.assertEqual(result, EXPECTED_RECORD_COUNT)

    def test_fetch_gb_records_ids(self):
        result = fetch_gb_records(LOCUS, TAXID, False)
        with EXPECT_RECORDS_JSON.open() as f:
            EXPECTED_RECORD_IDS = json.load(f)
        self.assertEqual(result, EXPECTED_RECORD_IDS)
