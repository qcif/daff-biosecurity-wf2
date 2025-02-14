import unittest
from fetch_gb_records import fetch_sources, fetch_gb_records

# Constants
EXPECTED_SOURCES = {'34577062': None}
EXPECTED_RECORD_IDS = ['2579155578']

gi_number = "34577062"
database = "nuccore"
ACCESSIONS = [gi_number]
LOCUS = 'COI'
TAXID = "9606"


class TestFetchRecords(unittest.TestCase):

    # @patch('fetch_gb_records.fetch_metadata')
    def test_fetch_sources(self):
        # pdb.set_trace()
        # mock_fetch_metadata.return_value = MOCK_METADATA

        result = fetch_sources(ACCESSIONS)
        print("ACTUAL SOURCES: ", result)
        print("EXPECTED SOURCES: ", EXPECTED_SOURCES)
        self.assertEqual(result, EXPECTED_SOURCES)

    # @patch('fetch_gb_records.fetch_')
    def test_fetch_gb_records(self):
        result = fetch_gb_records(LOCUS, TAXID, True)
        print("ACTUAL RECORD ID: ", result)
        print("EXPECTED RECORD ID: ", EXPECTED_RECORD_IDS)
        self.assertEqual(result, EXPECTED_RECORD_IDS)


if __name__ == '__main__':
    unittest.main()
