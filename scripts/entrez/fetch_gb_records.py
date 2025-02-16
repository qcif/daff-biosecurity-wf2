"""Use Entrez module from Biopython to fetch sequence data
from NCBI/Genbank."""

import sys
from Bio import Entrez
import re
from utils.config import Config
config = Config()

DEBUG_REQUESTS = True
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


def fetch_metadata(gis, database="nuccore"):
    metadata = []
    for i in range(0, len(gis), 10):
        batch = gis[i:i+10]
        ids = ",".join(batch)
        handle = Entrez.efetch(
            db=database,
            id=ids,
            rettype="gb",
            retmode="text")
        metadata_list = handle.read()
        handle.close()

        metadata.append(metadata_list)
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

    return publications


def fetch_sources(accessions) -> dict[str, dict]:
    '''Fetch a list of publication sources for each accession.'''
    sources = {}
    metadata_batches = fetch_metadata(accessions)
    for metadata in metadata_batches:
        formatted_metadata = metadata.split("//")
        for i, accession in enumerate(accessions):
            if i < len(formatted_metadata):
                publications = parse_metadata(formatted_metadata[i])
                sources[accession] = publications
    return sources


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
    payload_size = sys.getsizeof(results)
    print(f"[fetch_gb_records] Payload size: {payload_size} bytes")
    handle.close()
    if count:
        return results["Count"]
    return results["IdList"]


# Example usage: fetch FASTA and metadata -------------------------------------
# if __name__ == "__main__":
    # start_time = time.time()
    # gi_number = "34577062"
    # gi_number = "377581039"
    # accessions = [gi_number, '1066585321', '1393953329', '1519311736',
    #               '1519244926', '2462587281', '2462587279', '2462587277',
    #               '2462587276', '2462587274']
    # accessions = [gi_number, "377581039"]
    # database = "nuccore"

    # accession_metadata = fetch_metadata(gi_number)
    # publications_list = parse_metadata(accession_metadata)
    # print(f"[publications list]: {publications_list}")
    # payload_size = sys.getsizeof(accession_metadata)
    # print(f"[fetch_metadata] Payload size: {payload_size} bytes")
    # end_time = time.time()
    # execute_time = end_time - start_time
    # print(f"Request time: {execute_time:.2f} seconds")

    # source = fetch_sources(accessions)
    # print("[fetch sources]: ", source)

    # start_time = time.time()
    # taxid = 9606
    # locus = "COI"
    # record_count = fetch_gb_records(locus, taxid)
    # print("Record count: ", record_count)

    # end_time = time.time()
    # execute_time = end_time - start_time
    # print(f"Request time: {execute_time:.2f} seconds")
