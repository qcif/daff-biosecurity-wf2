from pathlib import Path
import requests
import logging
from xml.etree import ElementTree
from Bio import SeqIO
from src.utils import errors

logger = logging.getLogger(__name__)
ID_ENGINE_URL = "http://v4.boldsystems.org/index.php/Ids_xml"
BOLD_API_URL = "http://v4.boldsystems.org/index.php/API_Public/combined?"
WORK_CODE = 200
# query_dir = Path(__file__).resolve().parents[3] / 'output' / 'errors' / 'errors.json'


def main():
    fasta_file = (
        Path(__file__).resolve().parents[3]
        / 'tests'
        / 'test-data'
        / 'queries.fasta'
    )
    bold_taxa_engine = BoldTaxa(fasta_file)
    # print(f"Sequences extracted: {bold_taxa_engine.sequences}")

    # print(f"Hits Result: {bold_taxa_engine.hits_result}")

    # print(f"Extracted Taxa: {bold_taxa_engine.taxa}")

    print(f"Fetched Records: {bold_taxa_engine.records}")

    # taxon_counts = bold_taxa_engine.taxon_count()
    # taxon_collectors = bold_taxa_engine.taxon_collectors()
    # taxon_taxonomy = bold_taxa_engine.taxon_taxonomy()
    # print("Taxon Counts:", taxon_counts)
    # print("Taxon Collectors:", taxon_collectors)
    # print("Taxon Taxonomy:", taxon_taxonomy)


class BoldTaxa:
    """Fetch metadata for given taxa from the BOLD API."""
    def __init__(self, fasta_file: Path):
        self.fasta_file = fasta_file
        self.sequences = self._read_sequence_from_fasta(fasta_file)
        self.hits_result = self._bold_sequence_search(self.sequences)
        self.taxa = self._extract_taxa(self.hits_result)
        self.records = self._fetch_records(self.taxa)

    def taxon_count(self) -> dict[str, int]:
        """Return a count of each taxon occurences."""
        taxa_counts = {}
        for record in self.records:
            taxon = record.get("species_name", "")
            taxa_counts[taxon] = taxa_counts.get(taxon, 0) + 1
        return taxa_counts

    def taxon_collectors(self) -> dict[str, list[str]]:
        """Return a dictionary of taxa and related collectors."""
        taxa_collectors = {}
        for record in self.records:
            taxon = record.get("species_name", "")
            collector = record.get("collectors", "")
            if taxon:
                if taxon not in taxa_collectors:
                    taxa_collectors[taxon] = set()
                if collector:
                    taxa_collectors[taxon].update(collector.split(","))
        return taxa_collectors

    def taxon_taxonomy(self) -> dict[str, dict[str, str]]:
        """Build a taxonomy dictionary."""
        taxonomy_dict = {}
        for record in self.records:
            species_name = record.get("species_name", "")
            if species_name:
                taxonomy_dict[species_name] = {
                    "phylum": record.get("phylum_name", ""),
                    "class": record.get("class_name", ""),
                    "order": record.get("order_name", ""),
                    "family": record.get("family_name", ""),
                    "genus": record.get("genus_name", ""),
                    "species": species_name,
                }
        return taxonomy_dict

    def _read_sequence_from_fasta(self, fasta_file: Path):
        """Read sequence from fasta file."""
        sequences = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append(str(record.seq))
        return sequences

    def _bold_sequence_search(self, sequences: list[str],
                              db: str = "COX1_SPECIES_PUBLIC"
                              ) -> list[dict[str, any]]:
        """Submit a sequence search request to BOLD API."""
        hits_result = []
        for i, sequence in enumerate(sequences):
            params = {
                "sequence": sequence,
                "db": db
            }
            response = requests.get(ID_ENGINE_URL, params=params)
            if response.status_code >= 400:
                msg = (
                    f"Error for sequence {i + 1}: HTTP {response.status_code}"
                )
                logger.error(
                    f"Error for sequence {i + 1}: HTTP {response.status_code}"
                )
                errors.write(
                    errors.LOCATIONS.BOLD_ID_ENGINE,
                    msg,
                    None,
                    data={
                        "query_ix": i,
                        "response": response.text,
                    },
                    # query_dir=query_dir
                )
                continue

            root = ElementTree.fromstring(response.text)
            for match in root.findall("match"):
                result = {
                    "hit_id": match.find("ID").text,
                    "sequence_description": (
                        match.find("sequencedescription").text
                    ),
                    "database": match.find("database").text,
                    "citation": match.find("citation").text,
                    "taxonomic_identification": (
                        match.find("taxonomicidentification").text
                    ),
                    "similarity": match.find("similarity").text,
                    "specimen": match.find("specimen").text,
                    "url": match.find("specimen/url").text,
                    "country": (
                        match.find("specimen/collectionlocation/country").text
                    ),
                    "latitude": (
                        match.find(
                            "specimen/collectionlocation/coord/lat"
                        ).text
                        if match.find(
                            "specimen/collectionlocation/coord/lat"
                        ) is not None else ""
                    ),
                    "longitude": (
                        match.find(
                            "specimen/collectionlocation/coord/lon"
                        ).text
                        if match.find(
                            "specimen/collectionlocation/coord/lon"
                        ) is not None else ""
                    ),
                }
                hits_result.append(result)
        return hits_result

    def _extract_taxa(self, hits_result: list[dict]) -> list[str]:
        """Extract taxa (taxonomic_identification)."""
        taxa = []
        for hit in hits_result:
            taxonomic_identification = hit.get("taxonomic_identification")
            if (
                taxonomic_identification
                and taxonomic_identification not in taxa
            ):
                taxa.append(taxonomic_identification)
        return taxa

    def _fetch_records(self, taxa: list[str]) -> list[dict]:
        """Fetch accessions by calling BOLD public API with Taxa and
        save accessions to a .tsv file."""
        records = []
        taxa_param = "|".join(taxa)
        params = {
            "taxon": taxa_param,
            "format": "tsv"
        }

        response = requests.get(BOLD_API_URL, params=params)
        if response.status_code == WORK_CODE:
            lines = response.text.splitlines()
            headers = lines[0].split("\t")
            for line in lines[1:]:
                row = dict(zip(headers, line.split("\t")))
                records.append(row)
            logger.info(
                f"Records fetched successfully: {len(records)} records."
            )
        else:
            msg = (
                f"Error status code: {response.status_code}"
            )
            logger.error(f"Error status code: {response.status_code}")
            errors.write(
                errors.LOCATIONS.BOLD_TAXA,
                msg,
                None,
                data={
                    "taxa": taxa,
                    "response": response.text,
                },
                # query_dir=query_dir
            )

        return records


# Example usage
if __name__ == "__main__":
    main()
