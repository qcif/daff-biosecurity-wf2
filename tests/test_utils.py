import tempfile
import shutil
import unittest
from pathlib import Path

from src.utils.flags import Flag, FLAGS


class TestUtils(unittest.TestCase):

    TEST_SPECIES = 'Test species'
    TEST_TARGET_TYPE = 'candidate'

    def setUp(self):
        self.query_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.query_dir)

    def test_flags(self):
        Flag.write(self.query_dir, FLAGS.POSITIVE_ID, FLAGS.A)
        Flag.write(self.query_dir, FLAGS.TOI, FLAGS.B)
        Flag.write(
            self.query_dir,
            FLAGS.SOURCES,
            FLAGS.A,
            target=self.TEST_SPECIES,
            target_type=self.TEST_TARGET_TYPE,
        )
        flags = Flag.read(self.query_dir)
        flag_1 = flags[FLAGS.POSITIVE_ID]
        self.assertIsInstance(flag_1, Flag)
        self.assertEqual(flag_1.value, FLAGS.A)
        flag_2 = flags[FLAGS.TOI]
        self.assertIsInstance(flag_2, Flag)
        self.assertEqual(flag_2.value, FLAGS.B)
        flag_4 = flags[FLAGS.SOURCES][self.TEST_TARGET_TYPE][
            self.TEST_SPECIES]
        self.assertIsInstance(flag_4, Flag)
        self.assertEqual(flag_4.value, FLAGS.A)
        flags_json = Flag.read(self.query_dir, as_json=True)
        flag_4 = flags_json[FLAGS.SOURCES][self.TEST_TARGET_TYPE][
            self.TEST_SPECIES]
        self.assertIsInstance(flag_4, str)
        self.assertEqual(flag_4, FLAGS.A)
