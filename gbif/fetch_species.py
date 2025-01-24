"""Fetch species under the given taxon from GBIF API."""

import pygbif

LIMIT = 500


def fetch_species(taxon, rank, include_extinct=False):
    """Fetch species under the given taxon from GBIF API."""

    def filter_records(records):
        return [
            r for r in records
            if r['taxonomicStatus'] in ['ACCEPTED', 'DOUBTFUL']  #! correct?
        ]

    res = pygbif.species.name_suggest(
        q=taxon,
        rank=rank,
        limit=1,
    )
    genus_key = res[0]['genusKey']

    i = 0
    end_of_records = False
    records = []
    params = {
        'rank': 'species',
        'higherTaxonKey': genus_key,
        'limit': LIMIT,
        'offset': i * LIMIT,
    }
    if not include_extinct:
        params['isExtinct'] = False

    while not end_of_records:
        res = pygbif.species.name_lookup(**params)
        records += filter_records(res['results'])
        end_of_records = res['endOfRecords']
        i += 1

    if rank == 'species':
        records = [
            r for r in records
            if r['canonicalName'].lower() != taxon.lower()
        ]

    return records


def fetch_species_keys_for_country(
    taxon,
    rank,
    country_code,
):
    """Fetch species under the given taxon from GBIF API."""
    res = pygbif.species.name_suggest(
        q=taxon,
        rank=rank,
        limit=1,
    )
    genus_key = res[0]['genusKey']

    i = 0
    end_of_records = False
    records = []
    while not end_of_records:
        res = pygbif.occurrences.search(
            genusKey=genus_key,
            country=country_code,
            facet="speciesKey",
            facetLimit=LIMIT,
            offset=i * LIMIT,
            limit=1,  # don't need every occurence for each species
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
    records = []
    for species in species_counts:
        species_key = species.get("name")
        records.append(int(species_key))

    return records


if __name__ == '__main__':
    taxon = 'Phascolarctos cinereus'  # Koala
    country_code = 'AU'
    species_records = fetch_species(taxon, 'species')
    species_keys_for_country = fetch_species_keys_for_country(
        taxon, 'species', country_code)
    species_in_country = [
        species for species in species_records
        if species['speciesKey'] in species_keys_for_country
    ]
    print(f"Species found in country {country_code} related to {taxon}:")
    if not species_in_country:
        print("No related species found.")
    for species in species_in_country:
        print(species['canonicalName'])
