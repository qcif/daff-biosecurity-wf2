"""Fetch species under the given taxon from GBIF API."""

import argparse
import csv
import logging
import pygbif
import time
from pathlib import Path

LIMIT = 500
ACCEPTED_STATUS = ['ACCEPTED', 'DOUBTFUL']
GBIF_SPECIES_BASE_URL = 'https://www.gbif.org/species/'

logger = logging.getLogger(__name__)


def main():
    start_time = time.time()
    args = _parse_args()
    relatives = RelatedTaxaGBIF(args.taxon, args.rank)
    logger.info(f"Found {len(relatives.related_species)} species related to"
                f" {args.taxon}.")

    if args.country:
        species_in_country = relatives.for_country(args.country)
        logger.info(f"Species found in country {args.country} related to"
                    f" {args.taxon}:")
        if not species_in_country:
            logger.info("No related species found.")
        else:
            for species in species_in_country:
                logger.info(species['canonicalName'])

    _write_output(args.output, relatives.related_species, species_in_country)
    logger.info(f"Fetched {len(relatives.related_species)} species in"
                f" {time.time() - start_time:.2f} seconds")


def _parse_args():
    parser = argparse.ArgumentParser(
        description="Fetch species under the given taxon from GBIF API."
    )
    parser.add_argument(
        "taxon",
        type=str,
        help="Taxon to fetch species for."
    )
    parser.add_argument(
        "--rank",
        type=str,
        default="species",
        help="Rank of the taxon."
    )
    parser.add_argument(
        "--country",
        type=str,
        help="Country code to filter species by."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default="related.csv",
        help="File path to save related species to."
    )
    return parser.parse_args()


def _get_status(record):
    return record.get('status') or record.get('taxonomicStatus')


def _get_gbif_url(record):
    species_key = record.get('speciesKey')
    return f"{GBIF_SPECIES_BASE_URL}{species_key}"


def _write_output(
    filepath,
    related_species,
    species_in_country=None,
):
    fields = [
        'canonicalName',
        ('GBIF URL', _get_gbif_url),
        'kingdom',
        'phylum',
        'class',
        'order',
        'family',
        'genus',
        'species',
        ('status', _get_status),
        'rank',
        'synonym',
    ]
    logger.info(f"Writing related species to {filepath}...")
    with open(filepath, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for species in related_species:
            writer.writerow({
                k: species[k]
                if not isinstance(k, tuple)
                else k[1](species)
                for k in fields
            })

    if species_in_country:
        country_path = filepath.with_name(f"{filepath.stem}_country.csv")
        logger.info(f"Writing species for country to {country_path}...")
        with open(country_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            for species in species_in_country:
                writer.writerow({
                    k: species[k]
                    if not isinstance(k, tuple)
                    else k[1](species)
                    for k in fields
                })


class RelatedTaxaGBIF:
    """Fetch taxonomic relatives for a given taxon from GBIF API."""

    INCLUDE_EXTINCT = False

    def __init__(self, taxon, rank='species'):
        self.taxon = taxon
        self.rank = rank
        self.genus_key = self._get_genus_key(taxon, rank)
        self.related_species = self.fetch_related()

    def _get_genus_key(self, taxon, rank):
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
            record[status_key] in ACCEPTED_STATUS
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
            'limit': LIMIT,
            'offset': i * LIMIT,
        }

        while not end_of_records:
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
            res = pygbif.occurrences.search(
                genusKey=self.genus_key,
                country=country_code,
                facet="speciesKey",
                facetLimit=LIMIT,
                offset=i * LIMIT,
                limit=1,  # we don't need every occurence for each species
            )
            records += res['results']
            try:
                end_of_records = len(res['facets'][0]['counts']) < LIMIT
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


if __name__ == '__main__':
    main()

    # Tested:
    taxon = 'Phascolarctos cinereus'    # Koala (0 relatives)
    taxon = 'Froggattella latispina'    # Ant (1/1 relatives)
    taxon = 'Encoptarthria grisea'      # Spider (2/4 relatives)
    taxon = 'Cheiloxena aitori'         # Leaf beetle (7/7 relatives)
