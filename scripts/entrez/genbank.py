"""Use Entrez module from Biopython to fetch sequence data
from NCBI/Genbank."""

import logging
import time
from Bio import Entrez
import re
from filelock import FileLock

from ..utils.config import Config

config = Config()
filelock = FileLock(config.entrez_lock_file)
logger = logging.getLogger(__name__)

DEBUG_REQUESTS = True
LOCI = {  # TODO: update with DAFF
  'coi': ['COI', 'COX', 'CO1', 'Cytochrome oxidase subunit 1'],
  # others to be confirmed with DAFF
}

Entrez.email = config.USER_EMAIL


def fetch_entrez(
    endpoint=Entrez.efetch,
    db='nuccore',
    **kwargs,
) -> str:
    """Fetch data from NCBI Entrez database.
    Queue requests to avoid rate limits.
    Retry-on-failure.
    """
    def read(handle):
        if endpoint == Entrez.efetch:
            return handle.read()
        return Entrez.read(handle)

    handle = None
    retries = config.ENTREZ_MAX_RETRIES
    kwargs.update({
        "db": db,
    })
    while True:
        with filelock:
            # Ensure that requests are sent at max 10/sec
            time.sleep(0.1)
        try:
            handle = endpoint(**kwargs)
            data = read(handle)
            break
        except Exception as exc:
            logger.warning(f"Error fetching data from Entrez API:\n{exc}")
            if not retries:
                logger.error("Max retries reached. Exiting.")
                raise exc
            logger.info(f"Retrying {retries} more times.")
            retries -= 1
        finally:
            if handle:
                handle.close()
    return data


def fetch_fasta(identifier, **kwargs):
    return fetch_entrez(
        ids=identifier,
        rettype="fasta",
        retmode="text",
        **kwargs,
    )


def fetch_sources(accessions, **kwargs) -> list[dict]:
    accession_sources = {}
    for i in range(0, len(accessions), 10):
        batch = accessions[i:i+10]
        ids_str = ",".join(batch)
        metadata_batches = fetch_entrez(
            id=ids_str,
            rettype="gb",
            retmode="text",
            **kwargs,
        )
        metadata_list = [
            record.strip(' \n')
            for record in metadata_batches.split("\n//\n")
            if record.strip(' \n')
        ]
        for record in metadata_list:
            accession = match_accession_to_metadata(record, batch)
            sources = parse_metadata_sources(record)
            accession_sources[accession] = sources

    missing_accessions = set(accessions) - set(accession_sources.keys())
    if missing_accessions:
        logger.warning("[Fetch sources] no data was returned for the"
                       " following accessions. They cannot be included in"
                       " source diversity analysis:"
                       f" {missing_accessions}")
        for accession in missing_accessions:
            accession_sources[accession] = []

    return accession_sources


def match_accession_to_metadata(record, batch):
    """Find the accession that matched the given from metadata record."""
    matching_accession = [
        a for a in batch
        if a.lower() in record.lower()
    ]
    if matching_accession:
        return matching_accession[0]

    logger.warning("Accession number found in metadata record does not"
                   " match any accession number in the batch. This record can"
                   " not be used for source diversity analysis.")


def parse_metadata_sources(metadata):
    '''Extract authors, title, journal and year from data.'''
    publications = []
    current_publications = {}

    for line in metadata.split("\n"):
        if line.startswith("REFERENCE"):
            if current_publications:
                publications.append(current_publications)
            current_publications = {}
        elif line.startswith("  AUTHORS"):
            current_publications['authors'] = (
                line.replace("AUTHORS", "")
                .strip()
                .split(", ")
            )
        elif line.startswith("  TITLE"):
            current_publications['title'] = (
                line.replace("TITLE", "")
                .strip()
            )
        elif line.startswith("  JOURNAL"):
            journal_line = line.replace("JOURNAL", "").strip()
            current_publications['journal'] = journal_line
            match = re.search(r'\((\d{4})\)', journal_line)
            if match:
                current_publications['year'] = int(match.group(1))

    if current_publications:
        publications.append(current_publications)

    return publications


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
    results = fetch_entrez(
        endpoint=Entrez.esearch,
        term=query,
        retmax=max_results,
    )
    if count:
        return int(results["Count"])
    return results["IdList"]
