from pathlib import Path
import requests
from xml.etree import ElementTree
import csv


class BOLDIdEngine:
    ID_ENGINE_URL = "http://v4.boldsystems.org/index.php/Ids_xml"

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

    def write_hits_tsv(self, xml_data):
        """Parse XML response to extract hit data."""


# Example usage
if __name__ == "__main__":
    fasta_file = Path(__file__).parent / 'example_coi.fasta'
    idEngine = BOLDIdEngine()

    # Read sequence from FASTA file
    sequence = idEngine.read_sequence_from_fasta(fasta_file)
    print("Read Sequence:", sequence)

    response = idEngine.sequence_search(sequence)
    print("BOLD API Response:", response)
