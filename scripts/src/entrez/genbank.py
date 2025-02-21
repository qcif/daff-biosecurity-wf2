"""Use Entrez module from Biopython to fetch sequence data
from NCBI/Genbank."""

import logging
import time
from Bio import Entrez
from filelock import FileLock
from xml.etree import ElementTree as ET

from ..utils.config import Config

config = Config()
logger = logging.getLogger(__name__)

DEBUG_REQUESTS = True
EFETCH_BATCH_SIZE = 10
AUTOMATED_ANNOTATION_TAG = '##Genome-Annotation-Data-START##'
LOCI = {  # TODO: update with DAFF - maybe read from allowed_loci.txt file?
  'coi': ['COI', 'COX', 'CO1', 'Cytochrome oxidase subunit 1'],
  # others to be confirmed with DAFF
}

Entrez.email = config.USER_EMAIL


class GbRecordSource:
    """A source of sequence record.

    Can be either a publication, set of authors or automated annotation.
    Instances can be compared with the `shares_publication_with` method, which
    determines if two sources share a publication source. This is not done by
    exact matching of publication metadata, but by establishing a match
    between a single metadata field (in order of preference: title, authors,
    journal). If the source is automated (e.g. automated genome annotation),
    it will only match with other automated sources.
    """

    def __init__(self, accession: str):
        self.accession = accession
        self.is_automated = False
        self._publications = []

    def __hash__(self) -> int:
        return hash(str(self._publications))

    def __bool__(self) -> bool:
        return len(self._publications) > 0

    @property
    def publications(self) -> list[dict]:
        return [
            {
                'authors': authors,
                'title': title,
                'journal': journal,
            }
            for authors, title, journal in self._publications
        ]

    def to_json(self) -> dict:
        return {
            'is_automated': self.is_automated,
            'publications': self.publications,
        }

    def add_publication(
        self,
        authors: list[str] = None,
        title: str = None,
        journal: str = None,
    ):
        self._publications.append((
            authors,
            title,
            journal,
        ))

    def matches(self, other) -> bool:
        """Return True if other has a publication source that matches self."""
        if not isinstance(other, GbRecordSource):
            raise ValueError(
                "Cannot only comparte GbRecordSource with an object of the"
                f" same type. Received type '{type(other)}'")
        self_publications = self.get_publication_repr()
        other_publications = other.get_publication_repr()
        has_shared_source = self_publications == other_publications
        has_shared_publication = any(
            pub for pub in self_publications
            if pub in other_publications
        )
        if has_shared_publication or has_shared_source:
            return True
        return False

    def get_publication_repr(self) -> list[str]:
        """Return string representations of publications that can be used to
        assert fuzzy uniqueness against another (case-insensitive).
        """
        def get_repr(p):
            if authors := p.get('authors'):
                # Remove punctuation and spaces to allow for fuzzy matching
                # between records
                return ', '.join(
                    a.replace('.', '')
                    .replace(',', '')
                    .replace(' ', '')
                    .lower()
                    for a in authors
                )
            if title := p.get('title', ''):
                # Don't include "direct submission" titles - not useful
                if len(title) > 20 or 'direct submission' not in title.lower():
                    return title.lower()
            if journal := p.get('journal'):
                return journal.lower()

        if self.is_automated:
            return []
        return [
            get_repr(p)
            for p in self.publications
        ]


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
        with FileLock(config.entrez_lock_file):
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


def fetch_sources(accessions, **kwargs) -> dict[str, GbRecordSource]:
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
            accession_sources[accession] = GbRecordSource(accession)

    return accession_sources


def parse_metadata(xml_str) -> dict[str, GbRecordSource]:
    root = ET.fromstring(xml_str)
    records = {}
    for record in root.findall('.//GBSeq'):
        accession = getattr(
            record.find('GBSeq_primary-accession'), 'text', None)
        source = GbRecordSource(accession)
        for ref in record.findall('.//GBReference'):
            authors = [
                x.text for x in ref.findall('.//GBAuthor')
            ]
            title = getattr(ref.find('GBReference_title'), 'text', None)
            journal = getattr(ref.find('GBReference_journal'), 'text', None)
            source.add_publication(
                authors=authors,
                title=title,
                journal=journal,
            )
        if AUTOMATED_ANNOTATION_TAG in ET.tostring(record).decode('utf-8'):
            source.is_automated = True
        records[accession] = source
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
