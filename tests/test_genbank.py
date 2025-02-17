import json
import unittest
import threading
from queue import Queue
from pathlib import Path

from scripts.entrez.genbank import fetch_sources, fetch_gb_records

DATA_DIR = Path(__file__).parent / 'test-data'
EXPECTED_RECORD_COUNT = 1686
EXPECT_RECORDS_JSON = DATA_DIR / 'genbank_expect_ids.json'
EXPECT_SINGLE_SOURCE_JSON = DATA_DIR / 'genbank_single_source.json'
EXPECT_MULTIPLE_SOURCES_JSON = DATA_DIR / 'genbank_multiple_sources.json'
ACCESSIONS_LIST_FILE = DATA_DIR / 'accessions.txt'

ACCESSION_1 = "NM_001126"
ACCESSION_2 = "HQ621368"
DATABASE = "nuccore"
SINGLE_ACCESSION = [ACCESSION_1]
MULTIPLE_ACCESSIONS = [ACCESSION_1, ACCESSION_2]
LOCUS = 'COI'
TAXID = "9606"
THREAD_COUNT = 20  # Number of parallel requests


class TestFetchRecords(unittest.TestCase):

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

    def test_fetch_gb_records_count(self):
        result = fetch_gb_records(LOCUS, TAXID, True)
        self.assertEqual(result, EXPECTED_RECORD_COUNT)

    def test_fetch_gb_records_ids(self):
        result = fetch_gb_records(LOCUS, TAXID, False)
        with EXPECT_RECORDS_JSON.open() as f:
            EXPECTED_RECORD_IDS = json.load(f)
        self.assertEqual(result, EXPECTED_RECORD_IDS)

    def test_parallel_requests(self):
        """Test that parallel requests do not get blocked by the API."""
        def worker(accession):
            try:
                result = fetch_sources([accession])
                results_queue.put((accession, result))
            except Exception as e:
                results_queue.put((accession, str(e)))

        results_queue = Queue()
        threads = []
        accessions = [
            a.strip()
            for a in ACCESSIONS_LIST_FILE.read_text().splitlines()
            if a.strip()
        ]
        for acc in accessions[:THREAD_COUNT]:  # Use the first 20 accessions
            t = threading.Thread(target=worker, args=(acc,))
            t.start()
            threads.append(t)

        # Wait for all threads to complete
        for t in threads:
            t.join()

        # Collect results
        results = {}
        while not results_queue.empty():
            acc, result = results_queue.get()
            results.update(result)

        # Assert that all requests succeeded (i.e., no exceptions)
        for acc, result in results.items():
            self.assertIsInstance(result, list, f"Failed for {acc}: {result}")
        self.assertEqual(len(results), THREAD_COUNT)
