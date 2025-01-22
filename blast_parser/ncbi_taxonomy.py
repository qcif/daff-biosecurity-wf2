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

        # if result.returncode != 0:
        #     print(f"Error retrieving taxonomy details: {result.stderr}")
        #     return
        
        # print(f"Taxonkit result: {result.stdout}")
        
        taxon_details_list = []
        for line in result.stdout.strip().split('\n'):
            # print(f"Processing line: {line}")
            fields = line.split('\t')
            # print("FIELDS IS: ", len(fields))
            # print("FIELD [0]: ", fields[0])
            # print("FIELD [1]: ", fields[1])
            # print("FIELD [2]: ", fields[2])
            # print("FIELDS: ", fields)
            if len(fields) == 3:
                taxid, taxon_details, ranks = fields[0], fields[1], fields[2]
                # taxid, taxon_details, ranks = line.split('\t')
                lineage_list = taxon_details.split(';')
                ranks_list = ranks.split(';')
                taxonomy = {
                    rank: name for rank,
                    name in zip(ranks_list, lineage_list)
                    if rank in TAXONOMIC_RANKS
                }
                # print("TAXONOMY:", taxonomy)
                taxon_details_list.append((taxid, taxonomy))
            else:
                print(f"Warning: Unexpected format in line: {line}")
        # print(f"Extracted Taxon Details: {taxon_details_list}")
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
        # print(f"TaxIDs: {taxids}")
        taxonomy_details_list = cls._get_taxon_details(taxids)
        # print("TAXON DETAIL LIST: ", taxonomy_details_list)
        taxonomies = [
            cls(
                species=taxonomy.get('species'),
                taxonomy=taxonomy,
                taxid=taxid
            )
            for taxid, taxonomy in taxonomy_details_list
        ]

        taxonomies_as_dict = [tax.as_dict() for tax in taxonomies]
        # for tax in taxonomies:
        #     print(f"Species: {tax.species}, Taxonomy: {tax.taxonomy}, TaxID: {tax.taxid}")
        # print(f"TAXON FROM NCBI: {taxonomies_as_dict}")
        return taxonomies_as_dict
