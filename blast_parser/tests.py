import unittest
from unittest.mock import patch
from types import SimpleNamespace
if __name__ == '__main__':
    from ncbi_taxonomy import NCBITaxonomy
else:
    from .ncbi_taxonomy import NCBITaxonomy

BLASTDBCMD_STDOUT = '1529436\n2711157\n1529435\n'
TAXONKIT_STDOUT = """1529436	1529436	cellular organisms;Eukaryota;Opisthokonta;Metazoa;Eumetazoa;Bilateria;Deuterostomia;Echinodermata;Pelmatozoa;Crinoidea;Articulata;Comatulida;Comatulidae;Comatulinae;Anneissia;Anneissia japonica	no rank;superkingdom;clade;kingdom;clade;clade;clade;phylum;clade;class;subclass;order;family;subfamily;genus;species
2711157	2711157	cellular organisms;Eukaryota;Opisthokonta;Metazoa;Eumetazoa;Bilateria;Deuterostomia;Echinodermata;Pelmatozoa;Crinoidea;Articulata;Comatulida;Comatulidae;Comatulinae;Anneissia;Anneissia pinguis	no rank;superkingdom;clade;kingdom;clade;clade;clade;phylum;clade;class;subclass;order;family;subfamily;genus;species
1529435	1529435	cellular organisms;Eukaryota;Opisthokonta;Metazoa;Eumetazoa;Bilateria;Deuterostomia;Echinodermata;Pelmatozoa;Crinoidea;Articulata;Comatulida;Comatulidae;Comatulinae;Anneissia;Anneissia bennetti	no rank;superkingdom;clade;kingdom;clade;clade;clade;phylum;clade;class;subclass;order;family;subfamily;genus;species
"""


class TestNCBITaxonomy(unittest.TestCase):
    @patch('subprocess.run')
    def test_it_can_extract_taxonomic_data_for_accessions(self, mock_run):
        def mock_run_side_effect(args, **kwargs):
            if "blastdbcmd" in args:
                return SimpleNamespace(stdout=BLASTDBCMD_STDOUT, returncode=0)
            elif args[0] == "taxonkit":
                return SimpleNamespace(stdout=TAXONKIT_STDOUT, returncode=0)
            raise NotImplementedError(
                f"Command not implemented for mock: {args[0]}")

        # Assign the side effect function to the mock
        mock_run.side_effect = mock_run_side_effect

        accessions = ['A', 'B', 'C']
        db_name = 'input'
        expected_output = {
            'A': {
                "species": "Anneissia japonica",
                "taxonomy": {
                    "superkingdom": "Eukaryota",
                    "kingdom": "Metazoa",
                    "phylum": "Echinodermata",
                    "class": "Crinoidea",
                    "order": "Comatulida",
                    "family": "Comatulidae",
                    "genus": "Anneissia",
                    "species": "Anneissia japonica"
                },
                "taxid": "1529436"
            },
            'B': {
                "species": "Anneissia pinguis",
                "taxonomy": {
                    "superkingdom": "Eukaryota",
                    "kingdom": "Metazoa",
                    "phylum": "Echinodermata",
                    "class": "Crinoidea",
                    "order": "Comatulida",
                    "family": "Comatulidae",
                    "genus": "Anneissia",
                    "species": "Anneissia pinguis"
                },
                "taxid": "2711157"
            },
            'C': {
                "species": "Anneissia bennetti",
                "taxonomy": {
                    "superkingdom": "Eukaryota",
                    "kingdom": "Metazoa",
                    "phylum": "Echinodermata",
                    "class": "Crinoidea",
                    "order": "Comatulida",
                    "family": "Comatulidae",
                    "genus": "Anneissia",
                    "species": "Anneissia bennetti"
                },
                "taxid": "1529435"
            }
        }

        # import and run NCBITaxonomy.extract here
        output = NCBITaxonomy.extract(db_name, accessions)
        # Make some assertions on the return value here
        self.assertEqual(output, expected_output)


if __name__ == '__main__':
    unittest.main()
