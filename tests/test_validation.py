"""Test input validation functions."""

import logging
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from p0_validation import (
    TAXDB_EXPECT_FILES,
    _validate_fasta,
    _validate_metadata,
    _validate_metadata_sample_id,
    _validate_metadata_locus,
    _validate_metadata_preliminary_id,
    _validate_metadata_taxa_of_interest,
    _validate_metadata_country,
    _validate_metadata_host,
    _validate_taxdbs,
)
from src.utils.errors import FASTAFormatError, MetadataFormatError
from src.utils.config import Config

logging.disable(logging.CRITICAL)

TEST_DATA_DIR = Path(__file__).parent / 'test-data'
FASTA_VALID = TEST_DATA_DIR / 'queries.fasta'
FASTA_INVALID_RESIDUE = TEST_DATA_DIR / 'validation/invalid_residues.fasta'
FASTA_INVALID_LENGTH = TEST_DATA_DIR / 'validation/invalid_max_length.fasta'
FASTA_INVALID_COUNT = TEST_DATA_DIR / 'validation/invalid_max_count.fasta'
FASTA_INVALID_DUPLICATE_ID = (
    TEST_DATA_DIR / 'validation/invalid_duplicate.fasta')
METADATA_VALID = TEST_DATA_DIR / 'metadata.csv'
METADATA_INVALID_COLUMNS = (
    TEST_DATA_DIR / 'validation/metadata_invalid_columns.csv')
METADATA_MISSING_SAMPLE_ID = (
    TEST_DATA_DIR / 'validation/metadata_missing_sample_id.csv')
METADATA_MISSING_SEQUENCE_ID = (
    TEST_DATA_DIR / 'validation/metadata_missing_sequence_id.csv')
METADATA_EXPECT_IDS = [
    'LC438549.1',
    'ON075825.1',
    'PP466915.1',
    'JQ585746.1',
    'LC547004.1',
]
MOCK_TAXONKIT_DIR = Path(tempfile.gettempdir())
MOCK_TAXONKIT_DIR.mkdir(exist_ok=True, parents=True)


class ValidationTestCase(unittest.TestCase):

    def test_it_can_validate_fasta_input(self):
        _validate_fasta(FASTA_VALID)
        with self.assertRaises(FASTAFormatError):
            _validate_fasta(FASTA_INVALID_RESIDUE)
        with self.assertRaises(FASTAFormatError):
            _validate_fasta(FASTA_INVALID_LENGTH)
        with self.assertRaises(FASTAFormatError):
            _validate_fasta(FASTA_INVALID_COUNT)
        with self.assertRaises(FASTAFormatError):
            _validate_fasta(FASTA_INVALID_DUPLICATE_ID)

    def test_it_can_validate_metadata_csv(self):
        _validate_metadata(METADATA_VALID, METADATA_EXPECT_IDS)
        with self.assertRaises(MetadataFormatError):
            _validate_metadata(METADATA_INVALID_COLUMNS, METADATA_EXPECT_IDS)
        with self.assertRaises(MetadataFormatError):
            _validate_metadata(METADATA_MISSING_SAMPLE_ID, METADATA_EXPECT_IDS)
        with self.assertRaises(MetadataFormatError):
            _validate_metadata(METADATA_MISSING_SEQUENCE_ID,
                               METADATA_EXPECT_IDS)

    def test_it_can_validate_metadata_sample_id(self):
        _validate_metadata_sample_id('abcd_1234-blah.blah')
        with self.assertRaises(MetadataFormatError):
            _validate_metadata_sample_id('abcd_1234-blah blah blah')

    def test_it_can_validate_metadata_locus(self):
        _validate_metadata_locus('COI')
        os.environ['ALLOWED_LOCI_PATH'] = str(
            TEST_DATA_DIR / 'allowed_loci.txt')
        _validate_metadata_locus('coi')
        with self.assertRaises(MetadataFormatError):
            _validate_metadata_locus('GGG')

    def test_it_can_validate_metadata_preliminary_id(self):
        _validate_metadata_preliminary_id('Homo')
        _validate_metadata_preliminary_id('Homo sapiens')
        with self.assertRaises(MetadataFormatError):
            _validate_metadata_preliminary_id('Homo-sapiens')

    def test_it_can_validate_metadata_taxa_of_interest(self):
        _validate_metadata_taxa_of_interest('Homo sapiens')
        _validate_metadata_taxa_of_interest('Homo sapiens|Pan troglodytes')
        with self.assertRaises(MetadataFormatError):
            _validate_metadata_taxa_of_interest('Homo sapiens,Pan troglodytes')

    def test_it_can_validate_metadata_country(self):
        _validate_metadata_country('canada')
        with self.assertRaises(MetadataFormatError):
            _validate_metadata_country('canada and mexico')

    def test_it_can_validate_metadata_host(self):
        _validate_metadata_host('Cut flowers Rosa')

    @patch.object(Config, 'TAXONKIT_DATA', new=MOCK_TAXONKIT_DIR)
    def test_it_can_validate_taxdbs_path(self):
        for required_file in TAXDB_EXPECT_FILES:
            path = MOCK_TAXONKIT_DIR / required_file
            path.touch()
        _validate_taxdbs(MOCK_TAXONKIT_DIR)
        path.unlink()
        with self.assertRaises(FileNotFoundError):
            _validate_taxdbs(MOCK_TAXONKIT_DIR)
