#!/usr/bin/env python3

import importlib
import os
import shutil
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

TEMPDIR_PREFIX = "integration_test_"


def print_green(text: str):
    print(f"\033[32m{text}\033[0m")


class IntegrationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for var in (
            'USER_EMAIL',
            'NCBI_API_KEY',
            'TAXONKIT_DATA',
        ):
            if var not in os.environ:
                raise EnvironmentError(
                    f"Environment variable {var} is not set. "
                    "Please set it before running integration tests. You may"
                    " wish to set this in a venv/bin/activate script or in"
                    " your shell profile.")
        cls.project_root = Path(__file__).parents[2]
        cls.scripts_root = cls.project_root / "scripts"
        cls.python = cls.project_root / "venv" / "bin" / "python"
        cls.taxdump_dir = Path.home() / ".taxonkit"
        cls.test_case_root = (
            cls.project_root / "tests" / "test-data"
            / "integration")

    def setUp(self):
        """Clean up old temp dirs and create a new one."""
        if not os.getenv("KEEP_OUTPUTS") == "1":
            tmp_root = Path(tempfile.gettempdir())
            for old_wdir in tmp_root.glob(f"{TEMPDIR_PREFIX}*"):
                if old_wdir.is_dir():
                    shutil.rmtree(old_wdir, ignore_errors=True)
        self.wdir_root = Path(tempfile.mkdtemp(prefix=TEMPDIR_PREFIX))
        self.test_cases = []
        self.completed_tests = []

    def tearDown(self):
        """Check if all tests passed and clean up."""
        print("\nCompleted tests:")
        for test_case in self.completed_tests:
            print(f"  - {test_case.name}")

        if len(self.completed_tests) != len(self.test_cases):
            print(f"Test failed. Wdir has been retained: {self.wdir_root}")
            return

        print("Test passed.")
        if os.getenv("KEEP_OUTPUTS") == "1":
            print("KEEP_OUTPUTS=1; output dirs have been retained at:"
                  f" {self.wdir_root}")
        else:
            print(f"Cleaning up: {self.wdir_root}")
            shutil.rmtree(self.wdir_root, ignore_errors=True)

    def prepare_working_dir(self, test_case: Path) -> Path:
        """Copy test case files into a fresh working directory"""
        wdir = self.wdir_root / test_case.name
        shutil.copytree(test_case, wdir)
        return wdir

    def patch_and_run(self, module_name, patched_args):
        mock_args = Namespace(**patched_args)
        module = importlib.import_module(f"scripts.{module_name}")
        with patch.object(
            module,
            "_parse_args",
            return_value=mock_args,
        ):
            module.main()

    def test_integration_cases(self):
        test_cases = [
            path for path in sorted(self.test_case_root.iterdir())
            if path.is_dir()
        ]
        self.test_cases = test_cases
        for test_case in test_cases:
            limit_test_case = os.getenv("RUN_TEST_CASE")
            if (
                limit_test_case
                and limit_test_case != test_case.name
            ):
                print(
                    f"Skipping test case '{test_case.name}' - env var "
                    f" RUN_TEST_CASE={limit_test_case} has been set.")
                continue

            with self.subTest(test_case=test_case.name):
                query_dir = None
                wdir = self.prepare_working_dir(test_case)
                os.environ['INPUT_FASTA_FILEPATH'] = str(
                    wdir / "query.fasta")
                os.environ['INPUT_METADATA_CSV_FILEPATH'] = str(
                    wdir / "metadata.csv")

                self.patch_and_run(
                    "p0_validation",
                    {
                        "metadata_csv": wdir / "metadata.csv",
                        "query_fasta": wdir / "query.fasta",
                        "taxdb_dir": self.taxdump_dir,
                        "bold": False,
                    },
                )
                print_green(f"\nTest case {test_case.name}: P0 PASS\n")

                self.patch_and_run(
                    "p1_parse_blast",
                    {
                        "blast_xml_path": wdir / "blast_result.xml",
                        "output_dir": wdir,
                    },
                )
                print_green(f"\nTest case {test_case.name}: P1 PASS\n")

                self.patch_and_run(
                    "p2_extract_taxonomy",
                    {
                        "taxids_csv": wdir / "taxids.csv",
                        "output_dir": wdir,
                    },
                )
                print_green(f"\nTest case {test_case.name}: P2 PASS\n")

                query_dir = next(wdir.glob("query_001*"))

                self.patch_and_run(
                    "p3_assign_taxonomy",
                    {
                        "query_dir": query_dir,
                        "output_dir": wdir,
                        "bold": False,
                    },
                )
                print_green(f"\nTest case {test_case.name}: P3 PASS\n")

                candidates_count_file = next(
                    query_dir.glob("candidates_count.txt"))
                with open(candidates_count_file) as f:
                    candidates_count = int(f.read().strip())

                if candidates_count < 4:
                    self.patch_and_run(
                        "p4_source_diversity",
                        {
                            "query_dir": query_dir,
                            "output_dir": wdir,
                        },
                    )
                    print_green(f"\nTest case {test_case.name}: P4 PASS\n")

                    self.patch_and_run(
                        "p5_db_coverage",
                        {
                            "query_dir": query_dir,
                            "output_dir": wdir,
                            "bold": False,
                        },
                    )
                    print_green(f"\nTest case {test_case.name}: P4 PASS\n")

                else:
                    print_green(
                        f"\nTest case {test_case.name}: P4/P5 SKIPPED - "
                        f"{candidates_count=} > 3\n"
                    )

                # Copy newick tree into query dir
                nwk_file = next(wdir.glob("*.nwk"))
                shutil.copy2(nwk_file, query_dir)

                self.patch_and_run(
                    "p6_report",
                    {
                        "query_dir": query_dir,
                        "output_dir": wdir,
                        "bold": False,
                    },
                )
                print_green(f"\nTest case {test_case.name}: P6 PASS\n")

                self.completed_tests.append(test_case)


if __name__ == "__main__":
    unittest.main()
