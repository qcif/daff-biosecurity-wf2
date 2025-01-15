"""docstring"""
import subprocess

class NCBITaxonomy:
    """docstring"""

    def __init__(self, species: str, taxonomy: dict[str, str]):
        self.species = species
        self.taxonomy = taxonomy

    @staticmethod
    def _get_taxid_from_blastdbcmd(accessions: list[str],db_name: str) -> list[str]:
        '''Use blastdbcmd to retrieve the taxonomy ID associated with the accession number.'''
        accession_list = ",".join(accessions)
        result = subprocess.run(
            ['blastdbcmd', '-db', db_name, '-entry', accession_list, '-outfmt', '%T'],
            capture_output=True,
            text=True
        )

        taxids = result.stdout.strip().split('\n')
        return taxids


    def _get_lineage_from_taxonkit():
        '''Use taxonkit lineage to extract the taxonomy details.'''
        return 
    

    def _parse_taxon_details():
        '''The result written to stdout can be parsed in Python and taxonomy information extracted.'''
        return


    @classmethod
    def extract(cls, accessions: list[str]) -> list['NCBITaxonomy']:
        """Extract taxonomy information from NCBI given a list of accessions."""

        # Extract things
        db = "your_database_name"
        
        # Extract taxid
        taxids = cls._get_taxid_from_blastdbcmd(accessions, db)
        # Extract lineage

        # Parse lineage to taxonomy dict


        return [
            cls(
                species=None,
                taxonomy=None,
            )
            for accession in taxonomies  # TODO
        ]
    

    