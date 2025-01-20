"""docstring"""
import subprocess

TAXONOMIC_RANKS = [
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

    def __init__(self, species: str, taxonomy: dict[str, str], taxid: list[str]):
        self.species = species
        self.taxonomy = taxonomy
        self.taxid = taxid

    @staticmethod
    def _retrieve_taxid(accessions: list[str], db_name: str) -> list[str]:
        '''Use blastdbcmd to retrieve the taxonomy ID
        associated with the accession number.'''
        accession_list = ",".join(accessions)
        result = subprocess.run(
            ['blastdbcmd', '-db', db_name, '-entry',
                accession_list, '-outfmt', '%T'],
            capture_output=True,
            text=True
        )

        taxids = result.stdout.strip().split('\n')
        return taxids

    def _get_taxon_details(taxids: list[str]) -> list[dict[str, str]]:
        '''Use taxonkit lineage to extract the taxonomy details.'''
        taxid_list = "\n".join(taxids)
        result = subprocess.run(
            ['taxonkit', 'lineage', '-R'],
            input=taxid_list,
            capture_output=True,
            text=True
        )
        taxon_details = []
        for line in result.stdout.strip().split('\n'):
            taxid, taxon_details, ranks = line.split('  ')
            lineage_list = taxon_details.split(';')
            ranks_list = ranks.split(';')
            taxonomy = {
                rank: name for rank,
                name in zip(ranks_list, lineage_list)
                if rank in TAXONOMIC_RANKS
            }
            taxon_details.append((taxid, taxonomy))

        return taxon_details

    @classmethod
    def extract(cls, db, accessions: list[str]) -> list['NCBITaxonomy']:
        """Extract taxonomy information from NCBI
        given a list of accessions."""

        # Extract things
        taxids = cls._retrieve_taxid(accessions, db)
        taxonomy_details = cls._get_taxon_details(taxids)
        taxonomies = [
            cls(
                species=taxonomy.get('species'),
                taxonomy=taxonomy,
                taxid=taxid
            )
            for taxid, taxonomy in taxonomy_details
        ]

        return taxonomies
