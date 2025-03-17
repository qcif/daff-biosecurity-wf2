from pathlib import Path
import requests
from xml.etree import ElementTree
import csv


class BOLDIdEngineToAccessions:
    ID_ENGINE_URL = "http://v4.boldsystems.org/index.php/Ids_xml"
    BOLD_API_URL = "http://v4.boldsystems.org/index.php/API_Public/combined?"
    WORK_CODE = 200

    def read_sequence_from_fasta(self, filename):
        """Read sequence from fasta file."""
        with open(filename, "r") as file:
            lines = file.readlines()
            sequence = "".join(line.strip()
                               for line in lines
                               if not line.startswith(">"))
        return sequence

    def sequence_search(self, sequence, db="COX1_SPECIES_PUBLIC"):
        """Submit a sequence search request to BOLD API."""
        params = {
            "sequence": sequence,
            "db": db
        }
        response = requests.get(self.ID_ENGINE_URL, params=params)
        output_file = (Path(__file__).parent / 'output_response.xml')
        with open(output_file, "w") as file:
            file.write(response.text)
        return response.text

    def write_hits_tsv(self, xml_file: Path, hits_file_path: Path):
        """Parse XML response to extract hit data."""
        tree = ElementTree.parse(xml_file)
        root = tree.getroot()

        with open(hits_file_path, "w") as tsv_file:
            # Write the header
            tsv_file.write(
                "ID\tSequenceDescription\tDatabase\tCitation\t"
                "TaxonomicIdentification\tSimilarity\tCountry\tLatitude\t"
                "Longitude\tURL\n"
            )

            # Iterate over all matches
            for match in root.findall("match"):
                hit_id = match.find("ID").text
                sequence_description = match.find("sequencedescription").text
                database = match.find("database").text
                citation = match.find("citation").text
                taxonomic_identification = (
                    match.find("taxonomicidentification").text
                )
                similarity = match.find("similarity").text
                specimen = match.find("specimen")
                url = specimen.find("url").text
                collection_location = specimen.find("collectionlocation")
                country = collection_location.find("country").text
                coord = collection_location.find("coord")
                latitude = coord.find("lat").text
                longitude = coord.find("lon").text

                tsv_file.write(
                    f"{hit_id}\t{sequence_description}\t{database}\t"
                    f"{citation}\t{taxonomic_identification}\t"
                    f"{similarity}\t{country}\t{latitude}\t"
                    f"{longitude}\t{url}\n"
                )

        print(f"Records successfully written to {hits_file_path}")

    def extract_ids(self, hits_file_path: Path):
        """Extract IDs from hits.tsv."""
        ids = []
        with open(hits_file_path, "r") as file:
            reader = csv.DictReader(file, delimiter="\t")
            for row in reader:
                ids.append(row["ID"])
        return ids

    def fetch_ids_accessions(self, ids, ids_accessions_file_path: Path):
        """Fetch accessions by calling BOLD public API with IDs and
        save accessions to a .tsv file."""
        ids_param = "|".join(ids)
        params = {
            "ids": ids_param,
            "format": "tsv"
        }
        response = requests.get(self.BOLD_API_URL, params=params)
        if response.status_code == self.WORK_CODE:
            with open(ids_accessions_file_path, "w") as file:
                file.write(response.text)
                print(f"ID Accessions metadata saved to "
                      f"{ids_accessions_file_path}")
        else:
            print(f"Error status code: {response.status_code}")

    def extract_taxons(self, hits_file_path: Path):
        taxons = []
        with open(hits_file_path, "r") as file:
            reader = csv.DictReader(file, delimiter="\t")
            for row in reader:
                taxons.append(row["TaxonomicIdentification"])
        return taxons

    def fetch_taxons_accessions(self,
                                taxons,
                                taxons_accessions_file_path: Path):
        """Fetch accessions by calling BOLD public API with Taxons and
        save accessions to a .tsv file."""
        taxons_param = "|".join(taxons)
        params = {
            "taxon": taxons_param,
            "format": "tsv"
        }
        response = requests.get(self.BOLD_API_URL, params=params)
        if response.status_code == self.WORK_CODE:
            with open(taxons_accessions_file_path, "w") as file:
                file.write(response.text)
                print(f"Taxon Accessions metadata saved to"
                      f"{taxons_accessions_file_path}")
                # print("Accessions response:", response.text)
        else:
            print(f"Error status code: {response.status_code}")


# Example usage
if __name__ == "__main__":
    fasta_file = Path(__file__).parent / 'example_coi.fasta'
    id_engine_xml_path = Path(__file__).parent / "output_response.xml"
    idEngine = BOLDIdEngineToAccessions()

    # Read sequence from FASTA file
    sequence = idEngine.read_sequence_from_fasta(fasta_file)
    print("Read Sequence:", sequence)

    response = idEngine.sequence_search(sequence)
    print("BOLD API Response:", response)
    hits_file_path = Path(__file__).parent / "hits.tsv"
    idEngine.write_hits_tsv(id_engine_xml_path, hits_file_path)

    # extract ids to BOLD api for extracting metadata
    ids_accessions_file_path = Path(__file__).parent / "ids_accessions.tsv"
    taxons_accessions_path = Path(__file__).parent / "taxons_accessions.tsv"
    ids = idEngine.extract_ids(hits_file_path)
    taxons = idEngine.extract_taxons(hits_file_path)
    # print("ids: ", ids)

    idEngine.fetch_ids_accessions(ids, ids_accessions_file_path)
    idEngine.fetch_taxons_accessions(taxons, taxons_accessions_path)
