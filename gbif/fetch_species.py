"""Fetch species under the given taxon from GBIF API."""

import pygbif
import json
from time import time

LIMIT = 500


def fetch_species(taxon, rank):
    """Fetch species under the given taxon from GBIF API."""
    res = pygbif.species.name_lookup(
        q=taxon,
        rank=rank,
        limit=1,
    )
    genus_key = res['results'][0]['genusKey']

    i = 0
    end_of_records = False
    records = []
    while not end_of_records:
        res = pygbif.species.name_lookup(
            rank='species',
            higherTaxonKey=genus_key,
            limit=LIMIT,
            offset=i * LIMIT,
        )
        records += res['results']
        end_of_records = res['endOfRecords']
        i += 1

    if rank == 'species':
        records = [
            r for r in records
            if r['canonicalName'].lower() != taxon.lower()
        ]

    return records


if __name__ == '__main__':
    t0 = time()
    taxon = 'Poa annua'
    species = fetch_species(taxon, 'species')
    print(json.dumps(species[:3], indent=2))
    print(f"Found {len(species)} records in {time() - t0:.2f} seconds.")
