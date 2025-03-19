from pathlib import Path
import requests
from xml.etree import ElementTree
import csv

ID_ENGINE_URL = "http://v4.boldsystems.org/index.php/Ids_xml"
BOLD_API_URL = "http://v4.boldsystems.org/index.php/API_Public/combined?"
WORK_CODE = 200


def read_sequence_from_fasta(filename):
    """Read sequence from fasta file."""
    sequences = []
    current_sequence = []
    with open(filename, "r") as file:
        for line in file:
            if line.startswith(">"):
                if current_sequence:
                    sequences.append("".join(current_sequence))
                    current_sequence = []
            else:
                current_sequence.append(line.strip())
        if current_sequence:
            sequences.append("".join(current_sequence))
    return sequences


def save_sequences_to_tsv(sequences,
                          hits_file_path: Path,
                          db="COX1_SPECIES_PUBLIC"):
    """Submit a sequence search request to BOLD API and save reponse to tsv file."""
    with open(hits_file_path, "w") as tsv_file:
        tsv_file.write(
            "ID\tSequenceDescription\tDatabase\tCitation\t"
            "TaxonomicIdentification\tSimilarity\tCountry\tLatitude\t"
            "Longitude\tURL\n"
        )

        for i, sequence in enumerate(sequences):
            # Submit sequence to BOLD API
            params = {
                "sequence": sequence,
                "db": db
            }
            response = requests.get(ID_ENGINE_URL, params=params)
            if response.status_code != 200:
                print(f"Error for sequence {i + 1}: HTTP {response.status_code}")
                continue

            root = ElementTree.fromstring(response.text)
            for match in root.findall("match"):
                hit_id = match.find("ID").text
                sequence_description = match.find("sequencedescription").text
                database = match.find("database").text
                citation = match.find("citation").text
                taxonomic_identification = match.find("taxonomicidentification").text
                similarity = match.find("similarity").text
                specimen = match.find("specimen")
                url = specimen.find("url").text
                collection_location = specimen.find("collectionlocation")
                country = collection_location.find("country").text
                coord = collection_location.find("coord")
                latitude = coord.find("lat").text if coord is not None else ""
                longitude = coord.find("lon").text if coord is not None else ""

                # Write the match to the TSV file
                tsv_file.write(
                    f"{hit_id}\t{sequence_description}\t{database}\t"
                    f"{citation}\t{taxonomic_identification}\t"
                    f"{similarity}\t{country}\t{latitude}\t"
                    f"{longitude}\t{url}\n"
                )
    print(f"Results successfully written to {hits_file_path}")


def extract_ids(hits_file_path: Path):
    """Extract IDs from hits.tsv."""
    ids = []
    with open(hits_file_path, "r") as file:
        reader = csv.DictReader(file, delimiter="\t")
        for row in reader:
            ids.append(row["ID"])
    return ids


def fetch_ids_accessions(ids, ids_accessions_file_path: Path):
    """Fetch accessions by calling BOLD public API with IDs and
    save accessions to a .tsv file."""
    ids_counts = {}
    for id in ids:
        ids_counts[id] = 0
    ids_param = "|".join(ids)
    params = {
        "ids": ids_param,
        "format": "tsv"
    }
    response = requests.get(BOLD_API_URL, params=params)
    if response.status_code == WORK_CODE:
        with open(ids_accessions_file_path, "w") as file:
            file.write(response.text)
            print(f"ID Accessions metadata saved to "
                  f"{ids_accessions_file_path}")
        records = response.text.splitlines()

        # Count the number of times each taxon appears in the response
        for line in records[1:]:
            for id in ids:
                if id in line:
                    ids_counts[id] += 1

        print("Record count for each taxon:", records)
        for id, count in ids_counts.items():
            print(f"{id}: {count} records")
    else:
        print(f"Error status code: {response.status_code}")


def extract_taxons(hits_file_path: Path):
    taxons = []
    with open(hits_file_path, "r") as file:
        reader = csv.DictReader(file, delimiter="\t")
        for row in reader:
            taxons.append(row["TaxonomicIdentification"])
    return taxons


def fetch_taxons_accessions(taxons, taxons_accessions_file_path: Path):
    """Fetch accessions by calling BOLD public API with Taxons and
    save accessions to a .tsv file."""
    taxons_counts = {}
    for taxon in taxons:
        taxons_counts[taxon] = 0
    taxons_param = "|".join(taxons)
    params = {
        "taxon": taxons_param,
        "format": "tsv"
    }

    response = requests.get(BOLD_API_URL, params=params)
    if response.status_code == WORK_CODE:
        with open(taxons_accessions_file_path, "w") as file:
            file.write(response.text)
            print(f"Taxon Accessions metadata saved to"
                  f"{taxons_accessions_file_path}")
        records = response.text.splitlines()

        # Count the number of times each taxon appears in the response
        for line in records[1:]:  
            for taxon in taxons:
                if taxon in line: 
                    taxons_counts[taxon] += 1

        print("Record count for each taxon:", records)
        for taxon, count in taxons_counts.items():
            print(f"{taxon}: {count} records")

    else:
        print(f"Error status code: {response.status_code}")


# Example usage
if __name__ == "__main__":
    fasta_file = Path(__file__).resolve().parents[3] / 'tests' / 'test-data' / 'query.fasta'

    # Read sequence from FASTA file
    sequences = read_sequence_from_fasta(fasta_file)
    print("Read Sequence:", sequences)
    print("count of sequence: ", len(sequences))

    hits_file_path = Path(__file__).parent / "hits.tsv"
    save_sequences_to_tsv(sequences, hits_file_path)

    # extract ids to BOLD api for extracting metadata
    ids_accessions_file_path = Path(__file__).parent / "ids_accessions.tsv"
    taxons_accessions_path = Path(__file__).parent / "taxons_accessions.tsv"
    ids = extract_ids(hits_file_path)
    taxons = extract_taxons(hits_file_path)

    fetch_ids_accessions(ids, ids_accessions_file_path)
    fetch_taxons_accessions(taxons, taxons_accessions_path)
