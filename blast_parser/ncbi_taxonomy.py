"""docstring"""


class NCBITaxonomy:
    """docstring"""

    def __init__(self, species: str, taxonomy: dict[str, str]):
        self.species = species
        self.taxonomy = taxonomy

    @classmethod
    def extract(cls, accessions: list[str]) -> list['NCBITaxonomy']:
        """Extract taxonomy information from NCBI given a list of accessions."""

        # Extract things

        return [
            cls(
                species=None,
                taxonomy=None,
            )
            for accession in taxonomies  # TODO
        ]
