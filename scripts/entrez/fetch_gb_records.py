"""Example usage of the Entrez module from Biopython to fetch sequence data
from NCBI/Genbank."""

from Bio import Entrez
import re
from utils.config import Config
config = Config()

# Set your email (required by NCBI for API usage)
# This should be set from user input
# consider using a class to encapsulate this
Entrez.email = config.ENTREZ_EMAIL
# Entrez.email = "your_email@example.au"


# Fetch the FASTA sequence
def fetch_fasta(gi, database="nuccore"):
    handle = Entrez.efetch(db=database, id=gi, rettype="fasta", retmode="text")
    fasta_data = handle.read()
    handle.close()
    return fasta_data


# Fetch the accession metadata
def fetch_metadata(gi, database="nuccore"):
    handle = Entrez.efetch(db=database, id=gi, rettype="gb", retmode="text")
    metadata = handle.read()
    handle.close()
    return metadata


def search_nuccore_mrna(taxid, search_term):
    """Search the nuccore database for mRNA records."""
    query = (
        f'("{search_term}"[Gene name]'
        f' AND txid{taxid}[Organism])'
    )
    # Increase retmax for more results
    handle = Entrez.esearch(db="nuccore", term=query, retmax=100)
    search_results = Entrez.read(handle)
    handle.close()
    return search_results["IdList"]


def parse_metadata(metadata):
    '''Extract authors, title, journal and year from data.'''
    publications = []
    current_publications = {}

    for line in metadata.split("\n"):
        if line.startswith("REFERENCE"):
            if current_publications:
                publications.append(current_publications)
            current_publications = {}
        elif line.startswith("  AUTHORS"):
            current_publications['authors'] = line.replace("  AUTHORS   ", "")\
                .strip().split(", ")
        elif line.startswith("  TITLE"):
            current_publications['title'] = line.replace("  TITLE     ", "")\
                .strip()
        elif line.startswith("  JOURNAL"):
            journal_line = line.replace("  JOURNAL   ", "").strip()
            current_publications['journal'] = journal_line
            match = re.search(r'\((\d{4})\)', journal_line)
            if match:
                current_publications['year'] = int(match.group(1))

    if current_publications:
        publications.append(current_publications)

    return


def fetch_sources(accessions) -> dict[str, dict]:
    '''Fetch a list of publication srouces for each accession.'''
    sources = {}
    for accession in accessions:
        metadata = fetch_metadata(accession)
        publications = parse_metadata(metadata)
        sources[accession] = publications
    return sources


def fetch_gb_records(locus: str, taxid: int) -> list[dict]:
    '''Find matching GenBank records.'''
    search_term = locus
    record_ids = search_nuccore_mrna(taxid, search_term)
    return record_ids


# Example usage: fetch FASTA and metadata -------------------------------------

# GI number for the record you want to fetch
gi_number = "34577062"

# Database to search (e.g., "nuccore")
database = "nuccore"

# Fetch and display data
fasta_sequence = fetch_fasta(gi_number)
accession_metadata = fetch_metadata(gi_number)

print("FASTA Sequence:")
print(fasta_sequence)

print("\nAccession Metadata:")
print(accession_metadata)


# Example usage: keyword search -----------------------------------------------

# Taxonomic ID for the organism of interest
taxid = "9606"  # Replace with your desired taxid (e.g., for humans, 9606)

# Term to search (e.g., "COI")
search_term = "COI"

# Fetch IDs
record_ids = search_nuccore_mrna(taxid, search_term)

# Display results
if record_ids:
    print(f"Found {len(record_ids)} record(s):")
    print(record_ids)
else:
    print("No records found.")
