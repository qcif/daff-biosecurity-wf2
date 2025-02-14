"""Use Entrez module from Biopython to fetch sequence data
from NCBI/Genbank."""

import sys
import time
from Bio import Entrez
import re
from utils.config import Config
config = Config()

LOCI = {
  'coi': ['COI', 'COX', 'CO1', 'Cytochrome oxidase subunit 1'],
  # others to be confirmed with DAFF
}

Entrez.email = config.USER_EMAIL


def fetch_fasta(gi, database="nuccore"):
    handle = Entrez.efetch(db=database, id=gi, rettype="fasta", retmode="text")
    fasta_data = handle.read()
    handle.close()
    return fasta_data


def fetch_metadata(gi, database="nuccore"):
    handle = Entrez.efetch(db=database, id=gi, rettype="gb", retmode="text")
    metadata = handle.read()
    handle.close()
    return metadata


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
    '''Fetch a list of publication sources for each accession.'''
    sources = {}
    for accession in accessions:
        metadata = fetch_metadata(accession)
        publications = parse_metadata(metadata)
        sources[accession] = publications
    return sources


# TODO: NEED TO CHANGE AFTER DAFF MEETING
def fetch_gb_records(
    locus: str,
    taxid: int,
    count: bool = False,
) -> list[dict]:
    '''Find matching GenBank records.

    If count, returns a count of matching records, otherwise returns a list
    of matching accessions IDs.
    '''
    if locus.lower() in LOCI:
        gene_names = LOCI[locus.lower()]
    else:
        gene_names = [locus]
    query = ' OR '.join(
        [f'"{term}"[Gene name]' for term in gene_names])
    query += f' AND txid{taxid}[Organism])'
    max_results = 1 if count else 100
    handle = Entrez.esearch(db="nuccore", term=query, retmax=max_results)
    results = Entrez.read(handle)
    handle.close()
    if count:
        return results["Count"]
    return results["IdList"]


# Example usage: fetch FASTA and metadata -------------------------------------
if __name__ == "__main__":
    start_time = time.time()
    gi_number = "34577062"
    database = "nuccore"

    fasta_sequence = fetch_fasta(gi_number)
    accession_metadata = fetch_metadata(gi_number)

    print("FASTA Sequence:")
    print(fasta_sequence)

    print("\nAccession Metadata:")
    print(accession_metadata)

    # Example usage: keyword search Taxonomic ID for the organism of interest
    taxid = "9606"
    locus = "COI"
    record_ids = fetch_gb_records(locus, taxid, count=True)

    payload_size = sys.getsizeof(record_ids)
    print(f"Payload size: {payload_size} bytes")

    # Display results
    if record_ids:
        print(record_ids)
    else:
        print("No records found.")

    end_time = time.time()
    execute_time = end_time - start_time
    print(f"Request time: {execute_time:.2f} seconds")
