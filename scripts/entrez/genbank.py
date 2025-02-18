"""Use Entrez module from Biopython to fetch sequence data
from NCBI/Genbank."""

import logging
import time
from Bio import Entrez
from filelock import FileLock
from xml.etree import ElementTree as ET

try:
    from ..utils.config import Config
except ImportError:
    from utils.config import Config

config = Config()
filelock = FileLock(config.entrez_lock_file)
logger = logging.getLogger(__name__)

DEBUG_REQUESTS = True
EFETCH_BATCH_SIZE = 10
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
            time.sleep(0.11)
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
    for i in range(0, len(accessions), EFETCH_BATCH_SIZE):
        batch = accessions[i:i + EFETCH_BATCH_SIZE]
        ids_str = ",".join(batch)
        metadata_xml = fetch_entrez(
            id=ids_str,
            rettype="gb",
            retmode="xml",
            **kwargs,
        ).decode('utf-8')
        sources = parse_metadata(metadata_xml)
        accession_sources.update(sources)

    missing_accessions = set(accessions) - set(accession_sources.keys())
    if missing_accessions:
        logger.warning("[Fetch sources] no data was returned for the"
                       " following accessions. They cannot be included in"
                       " source diversity analysis:"
                       f" {missing_accessions}")
        for accession in missing_accessions:
            accession_sources[accession] = []

    return accession_sources


def parse_metadata(xml_str):
    root = ET.fromstring(xml_str)
    records = {}
    for record in root.findall('.//GBSeq'):
        accession = getattr(
            record.find('GBSeq_primary-accession'), 'text', None)
        publications = []
        for ref in record.findall('.//GBReference'):
            authors = [
                x.text for x in ref.findall('.//GBAuthor')
            ]
            title = getattr(ref.find('GBReference_title'), 'text', None)
            journal = getattr(ref.find('GBReference_journal'), 'text', None)
            publications.append({
                'authors': authors,
                'title': title,
                'journal': journal,
            })
        records[accession] = publications
    return records


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


if __name__ == '__main__':
    from pathlib import Path
    accessions = Path(
        '../tests/test-data/accessions.txt').read_text().splitlines()
    sources = fetch_sources(accessions[:20])
    print(f"Fetched sources for {len(sources)} accessions.")
