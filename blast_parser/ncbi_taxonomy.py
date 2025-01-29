import os
import subprocess
import tempfile

TAXONOMIC_RANKS = [
    'superkingdom',
    'kingdom',
    'phylum',
    'class',
    'order',
    'family',
    'genus',
    'species',
]


class NCBITaxonomy:
    """Use retrieving taxonomy ID to extract taxonomy details."""

    def __init__(
            self,
            species: str,
            taxonomy: dict[str, str],
            taxid: list[str]
    ):
        self.species = species
        self.taxonomy = taxonomy
        self.taxid = taxid

    @staticmethod
    def _retrieve_taxid(accessions: list[str], db_name: str) -> list[str]:
        '''Use blastdbcmd to retrieve the taxonomy ID
        associated with the accession number.'''
        accession_list = ",".join(accessions)
        result = subprocess.run(
            [
                'blastdbcmd',
                '-db', db_name,
                '-entry', accession_list,
                '-outfmt', '%T',
            ],
            capture_output=True,
            text=True
        )

        taxids = result.stdout.strip().split('\n')
        return taxids

    @staticmethod
    def _get_taxon_details(taxids: list[str]) -> list[dict[str, str]]:
        '''Use taxonkit lineage to extract the taxonomy details.'''

        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_file:
            temp_file.write("\n".join(taxids))
            temp_file.flush()
            temp_file_name = temp_file.name
            result = subprocess.run(
                [
                    'taxonkit',
                    'lineage',
                    '-R',
                    '-c', temp_file_name,
                ],
                capture_output=True,
                text=True
            )
        os.remove(temp_file_name)

        taxon_details_list = []
        for line in result.stdout.strip().split('\n'):
            fields = line.split('\t')[1:]
            if len(fields) == 3:
                taxid, taxon_details, ranks = fields[0], fields[1], fields[2]
                lineage_list = taxon_details.split(';')
                ranks_list = ranks.split(';')
                taxonomy = {
                    rank: name for rank,
                    name in zip(ranks_list, lineage_list)
                    if rank in TAXONOMIC_RANKS
                }
                taxon_details_list.append((taxid, taxonomy))
            else:
                print(f"Warning: Unexpected format in line: {line}")
        return taxon_details_list

    def as_dict(self) -> dict:
        """Convert The NCBITaxonomy to dictionary format."""
        return {
            "species": self.species,
            "taxonomy": self.taxonomy,
            "taxid": self.taxid
        }

    @classmethod
    def extract(cls, db, accessions: list[str]) -> dict[str, dict]:
        """Extract taxonomy information from NCBI
        given a list of accessions."""

        # Extract things
        taxids = cls._retrieve_taxid(accessions, db)
        taxonomy_details_list = cls._get_taxon_details(taxids)
        taxonomies = {
            accession: cls(
                species=taxonomy.get('species'),
                taxonomy=taxonomy,
                taxid=taxid
            ).as_dict()
            for accession, (taxid, taxonomy) in zip(
                accessions, taxonomy_details_list
            )
        }
        return taxonomies
