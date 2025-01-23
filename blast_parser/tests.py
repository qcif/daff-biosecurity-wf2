import unittest
from unittest.mock import patch
from types import SimpleNamespace

from ncbi_taxonomy import NCBITaxonomy

BLASTDBCMD_STDOUT = '?'
TAXONKIT_STDOUT = '?'


class TestSubprocessRun(unittest.TestCase):
    @patch('subprocess.run')
    def test_subprocess_run(self, mock_run):
        def mock_run_side_effect(args, **kwargs):
            if args[0] == "blastdbcmd":
                return SimpleNamespace(stdout=BLASTDBCMD_STDOUT, returncode=0)
            elif args[0] == "taxonkit":
                return SimpleNamespace(stdout=TAXONKIT_STDOUT, returncode=0)
            raise NotImplementedError(
                f"Command not implemented for mock: {args[0]}")

        # Assign the side effect function to the mock
        mock_run.side_effect = mock_run_side_effect

        # import and run NCBITaxonomy.extract here

        # Make some assertions on the return value here
