"""Fetch species under the given taxon from GBIF API."""

import logging
import filelock
import pygbif
import time

from ..utils import config

GBIF_SPECIES_BASE_URL = 'https://www.gbif.org/species/'

logger = logging.getLogger(__name__)
config = config.Config()


class Endpoint:
    """Define configuration for requesting different GBIF endpoints."""

    FAST = {  # occurence and name_suggest - 12 req/sec
        'interval_sec': 0.1,
        'lock_file': config.gbif_fast_lock_file,
    }
    SLOW = {  # name_lookup - 1.25 req/sec
        'interval_sec': 0.8,
        'lock_file': config.gbif_slow_lock_file,
    }

    def __init__(self, endpoint: dict):
        self.interval_sec = endpoint['interval_sec']
        self.lock_file = endpoint['lock_file']


class Throttle:
    """Use filelock to throttle requests to the GBIF API.

    This is necessary to avoid hitting the rate limit of the GBIF API, or
    overwhelming the server.

    A custom interval and lock file can be set to allow for different
    throttling intervals for each endpoint.
    """

    def __init__(
        self,
        endpoint: Endpoint,
    ):
        self.interval_sec = endpoint.interval_sec
        self.lock = filelock.FileLock(endpoint.lock_file)

    def __enter__(self):
        self.lock.acquire()
        time.sleep(self.interval_sec)
        self.lock.release()

    def __exit__(self, exc_type, exc_value, traceback):
        pass


class RelatedTaxaGBIF:
    """Fetch taxonomic relatives for a given taxon from GBIF API."""

    INCLUDE_EXTINCT = False

    def __init__(self, taxon, rank='species'):
        self.taxon = taxon
        self.rank = rank
        self.genus_key = self._get_genus_key(taxon, rank)
        self.related_species = self.fetch_related()

    def _get_genus_key(self, taxon, rank):
        endpoint = Endpoint(Endpoint.FAST)
        throttle = Throttle(endpoint)
        with throttle:
            res = pygbif.species.name_suggest(
                q=taxon,
                rank=rank,
                limit=20,
            )
        for record in res:
            if self._is_accepted(record) and 'genusKey' in record:
                logger.info(f"Genus key found for taxon {taxon}:"
                            f" {record['genusKey']}")
                return record['genusKey']
        raise ValueError(f"Genus key not found for taxon {taxon}.")

    def _is_accepted(self, record):
        status_key = 'status' if 'status' in record else 'taxonomicStatus'
        return (
            record[status_key] in config.GBIF_ACCEPTED_STATUS
            and (self.INCLUDE_EXTINCT or record.get('isExtinct') is not True)
        )

    def _filter_records(self, records):
        return [r for r in records if self._is_accepted(r)]

    def fetch_related(self):
        i = 0
        end_of_records = False
        records = []
        params = {
            'rank': 'species',
            'higherTaxonKey': self.genus_key,
            'limit': config.GBIF_LIMIT_RECORDS,
            'offset': i * config.GBIF_LIMIT_RECORDS,
        }

        while not end_of_records:
            endpoint = Endpoint(Endpoint.SLOW)
            throttle = Throttle(endpoint)
            with throttle:
                res = pygbif.species.name_lookup(**params)
                records += self._filter_records(res['results'])
                end_of_records = res['endOfRecords']
                i += 1

        if self.rank == 'species':
            records = [
                r for r in records
                if r['canonicalName'].lower() != self.taxon.lower()
            ]

        return records

    def for_country(self, country_code):
        i = 0
        end_of_records = False
        records = []
        while not end_of_records:
            endpoint = Endpoint(Endpoint.SLOW)
            throttle = Throttle(endpoint)
            with throttle:
                res = pygbif.occurrences.search(
                    genusKey=self.genus_key,
                    country=country_code,
                    facet="speciesKey",
                    facetLimit=config.GBIF_LIMIT_RECORDS,
                    offset=i * config.GBIF_LIMIT_RECORDS,
                    limit=1,  # don't need every occurence for each species
                )
            records += res['results']
            try:
                end_of_records = (
                    len(res['facets'][0]['counts'])
                    < config.GBIF_LIMIT_RECORDS)
            except (KeyError, IndexError):
                end_of_records = True
            i += 1

        species_facets = res.get("facets", [])
        species_counts = (
            species_facets[0].get("counts", [])
            if species_facets
            else []
        )

        # Retrieve species names for unique speciesKeys
        species_keys = []
        for species in species_counts:
            species_key = species.get("name")
            if species_key:
                species_keys.append(int(species_key))

        return [
            r for r in self.related_species
            if r['speciesKey'] in species_keys
        ]
