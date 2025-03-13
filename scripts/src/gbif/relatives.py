"""Fetch species under the given taxon from GBIF API."""

import logging
import filelock
import pygbif
import time
from pprint import pformat
from types import SimpleNamespace
from functools import cached_property

from ..utils import config
from ..utils.errors import APIError

logger = logging.getLogger(__name__)
config = config.Config()

MODULE_NAME = 'GBIF API'
GBIF_SPECIES_BASE_URL = 'https://www.gbif.org/species/'


class GBIFRecordNotFound(Exception):
    pass


class ENDPOINTS:
    """Define configuration for throttling different GBIF endpoints."""

    def __init__(self):
        self.FAST = SimpleNamespace(  # occurence and name_suggest - 12 req/sec
            interval_sec=0.1,
            lock_file=config.gbif_fast_lock_file,
        )
        self.SLOW = SimpleNamespace(  # name_lookup - 1.25 req/sec
            interval_sec=0.8,
            lock_file=config.gbif_slow_lock_file,
        )


class Throttle:
    """Use filelock to throttle requests to the GBIF API.

    This is necessary to avoid hitting the rate limit of the GBIF API, or
    overwhelming the server.

    A custom interval and lock file can be set to allow for different
    throttling intervals for each endpoint.
    """

    def __init__(
        self,
        endpoint: SimpleNamespace,  # one of ENDPOINTS
    ):
        self.interval_sec = endpoint.interval_sec
        self.lock = filelock.FileLock(endpoint.lock_file)

    def __enter__(self):
        self.lock.acquire()
        time.sleep(self.interval_sec)
        self.lock.release()

    def __exit__(self, exc_type, exc_value, traceback):
        pass

    def with_retry(self, func, args):
        retries = config.GBIF_MAX_RETRIES
        while True:
            try:
                with self:
                    return func(**args)
            except Exception as exc:
                retries -= 1
                if retries <= 0:
                    raise APIError(
                        'Failed to fetch data from GBIF API after'
                        f' {config.GBIF_MAX_RETRIES} retries. Please try'
                        f' resuming this job at a later time.'
                        f'\nException: {exc}'
                    )
                logger.warning(
                    "Exception encountered in call to GBIF endpoint"
                    f" {func.__name__} Retrying {retries} more times."
                    f" Exception: {exc}\n"
                    f" Args:\n{pformat(args)}")
                time.sleep(1)


class RANK:
    SPECIES = 1
    GENUS = 2
    FAMILY = 3
    ORDER = 4
    CLASS = 5
    PHYLUM = 6
    KINGDOM = 7

    @classmethod
    def from_string(cls, rank: str) -> str:
        return getattr(cls, rank.upper(), None)


class RelatedTaxaGBIF:
    """Fetch taxonomic relatives for a given taxon from GBIF API."""

    INCLUDE_EXTINCT = False

    def __init__(self, taxon):
        self.taxon = taxon
        self.record = self._get_taxon_record(taxon)
        self.key = self.record.get('key')
        self.genus_key = self.record.get('genusKey')
        self.rank = RANK.from_string(self.record.get('rank'))

    def __str__(self):
        return f"{self.__class__.__name__}: {self.taxon} ({self.rank})"

    def __repr__(self):
        return self.__str__()

    def _get_taxon_record(self, taxon):
        args = {
            'q': taxon,
            'limit': 20,
        }
        endpoints = ENDPOINTS()
        throttle = Throttle(endpoints.FAST)
        res = throttle.with_retry(
            pygbif.species.name_suggest,
            args,
        )
        for record in res:
            if self._is_accepted(record):
                logger.info(f"[{MODULE_NAME}] Record found for taxon"
                            f" '{taxon}' - rank:{record['rank']}"
                            f" genusKey:{record.get('genusKey')}")
                return record
        raise GBIFRecordNotFound(
            f"No GBIF record found for '{taxon}'. Taxonomic records cannot"
            " be retrieved for this species name - please check that this"
            " species name is correct.")

    def _is_accepted(self, record):
        status_key = 'status' if 'status' in record else 'taxonomicStatus'
        return (
            record[status_key] in config.GBIF_ACCEPTED_STATUS
            and (self.INCLUDE_EXTINCT or record.get('isExtinct') is not True)
        )

    def _filter_records(self, records):
        return [
            r for r in records
            if self._is_accepted(r)
            and 'canonicalName' in r
        ]

    @cached_property
    def relatives(self):
        """Fetch related species with self.genus_key."""
        i = 0
        end_of_records = False
        records = []
        args = {
            'rank': 'species',
            'higherTaxonKey': self.genus_key,
            'limit': config.GBIF_LIMIT_RECORDS,
        }

        while not end_of_records:
            args['offset'] = i * config.GBIF_LIMIT_RECORDS
            endpoints = ENDPOINTS()
            throttle = Throttle(endpoints.SLOW)
            res = throttle.with_retry(
                pygbif.species.name_lookup,
                args,
            )
            if i > 5:
                first_name = self._filter_records(
                    res['results'])[0]['canonicalName']
                if [
                    r for r in records
                    if r['canonicalName'] == first_name
                ]:
                    logger.warning(
                        'GBIF API claims endofRecords=False after >5 requests,'
                        ' but an existing canonicalName was returned.'
                        ' Exiting early to avoid infinite loop '
                        f' - {len(records)} records have been fetched.'
                        f' Taxon: {self.taxon}, Genus key: {self.genus_key}.'
                    )
            records += self._filter_records(res['results'])
            end_of_records = res['endOfRecords']
            i += 1

        if self.rank == RANK.SPECIES:
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
            args = {
                'genusKey': self.genus_key,
                'country': country_code,
                'facet': "speciesKey",
                'facetLimit': config.GBIF_LIMIT_RECORDS,
                'offset': i * config.GBIF_LIMIT_RECORDS,
                'limit': 1,  # don't need every occurence for each species
            }
            endpoints = ENDPOINTS()
            throttle = Throttle(endpoints.FAST)
            res = throttle.with_retry(
                pygbif.occurrences.search,
                args,
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
            r for r in self.relatives
            if r['speciesKey'] in species_keys
        ]
