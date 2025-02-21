"""Analyze the database coverage of target species at the given locus.

Database coverage is analysed at three levels:

1. Target species coverage: The number of records for the target species
2. Related species coverage: The number of records for species related to the
   target species
3. Related species from sample country of origin: as for (2), but only for
species which have occurence records in the same country as the target species.
"""

import argparse
import json
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed

from src.entrez import genbank
from src.utils import errors
from src.gbif.relatives import RelatedTaxaGBIF
from src.taxonomy import extract
from src.utils import existing_path, serialize
from src.utils.config import Config

logger = logging.getLogger(__name__)
config = Config()

MODULE_NAME = "Database Coverage"


def main():
    args = parse_args()
    config.configure(args.output_dir)
    config.configure_query_logger(args.query_dir)
    assess_coverage(args.query_dir)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "query_dir", type=existing_path, help="Path to query output directory")
    parser.add_argument(
        "--output_dir",
        type=existing_path,
        default=config.output_dir,
        help=f"Path to output directory. Defaults to {config.output_dir}.")
    return parser.parse_args()


def read_candidate_species(query_dir):
    candidates = config.read_json(query_dir / config.CANDIDATES_JSON)
    return candidates["species"]


def assess_coverage(query_dir):
    def get_args(func, query_dir, target, taxid, locus, country):
        if func == get_target_coverage:
            return func, taxid, locus
        elif func == get_related_coverage:
            return func, target, locus, query_dir
        elif func == get_related_country_coverage:
            return func, target, locus, country, query_dir

    locus = config.get_locus_for_query(query_dir)
    country = config.get_country_for_query(query_dir)
    candidates = read_candidate_species(query_dir)
    pmi = config.get_pmi_for_query(query_dir)
    toi_list = config.read_taxa_of_interest(query_dir)
    targets = candidates + toi_list + [pmi]
    target_taxids = extract.taxids(targets)
    taxid_to_species = {v: k for k, v in target_taxids.items()}

    target_gbif_taxa = {
        target: RelatedTaxaGBIF(target)
        for target in targets
    }

    tasks = [
        get_args(
            func,
            query_dir,
            target_gbif_taxa[target],
            taxid,
            locus,
            country,
        )
        for target, taxid in target_taxids.values()
        for func in (
            get_target_coverage,
            get_related_coverage,
            get_related_country_coverage,
        )
    ]

    with ThreadPoolExecutor() as executor:
        future_to_task = {
            executor.submit(*task): task
            for task in tasks
        }

    results = {}
    for future in as_completed(future_to_task):
        func, target = future_to_task[future][:2]
        try:
            if func.__name__ not in results:
                results[func.__name__] = {}
            results[func.__name__][target] = future.result()
        except Exception as exc:
            logger.error(
                f"[{MODULE_NAME}]: Error processing {func.__name__} for"
                f" {target}:\n{exc}")
            species_name = (
                taxid_to_species[target]
                if target.isdigit()
                else target
            )
            target_source = (
                "candidate" if species_name in candidates
                else "taxon of interest" if species_name in toi_list
                else "preliminary ID"
            )
            msg = (
                f"Error processing {func.__name__} for target species"
                f" '{species_name}' ({target_source}). This target could not"
                f" be evaluated.")
            errors.write(
                errors.LOCATIONS.DATABASE_COVERAGE,
                msg,
                exc,
                query_dir=query_dir)

    candidate_results = {}
    toi_results = {}
    pmi_results = {}

    for taxid in target_taxids.values():
        species_name = taxid_to_species[taxid]
        result = results[get_target_coverage.__name__][taxid]
        if species_name in candidates:
            candidate_results[species_name] = result
        if species_name in toi_list:
            toi_results[species_name] = result
        if species_name == pmi:
            pmi_results[species_name] = result

    write_db_coverage(
        query_dir,
        candidate_results,
        toi_results,
        pmi_results,
    )


def get_target_coverage(taxid, locus):
    """Return a count of the number of accessions for the given target."""
    # TODO: potential for caching taxid result here
    return genbank.fetch_gb_records(taxid, locus, count=True)


def get_related_coverage(gbif_target, locus, query_dir):
    """Return a count of the number of related species (same genus) and the
    number of species which have at least one accession in the database.
    """
    species_names = [
        r["canonicalName"]
        for r in gbif_target.related_species
    ]
    # TODO: potential for caching GBIF related taxa here
    results, errors = fetch_gb_records_for_species(species_names, locus)
    if errors:
        for species, exc in errors.items():
            msg = (
                f"Error fetching related species records from Entrez API:\n"
                f"(species: '{species}').")
            errors.write(
                errors.DB_COVERAGE_RELATED,
                msg,
                exc,
                query_dir=query_dir)
    return results


def get_related_country_coverage(
    gbif_target,
    locus,
    country,
    query_dir,
):
    species_names = [
        r["canonicalName"]
        for r in gbif_target.for_country(country)
    ]
    # TODO: potential for caching GBIF related/country taxa here
    results, errors = fetch_gb_records_for_species(species_names, locus)
    if errors:
        for species, exc in errors.items():
            msg = (
                f"Error fetching related/country species records from Entrez"
                f" API (species: '{species}').")
            errors.write(
                errors.DB_COVERAGE_RELATED_COUNTRY,
                msg,
                exc,
                query_dir=query_dir)
    return results


def fetch_gb_records_for_species(species_names, locus):
    """Fetch a count of the number of Genbank accessions for each species in
    the list.
    """
    taxids = extract.taxids(species_names)
    taxid_to_species = {
        v: k
        for k, v in taxids.items()
    }
    tasks = [
        (taxid, locus)
        for taxid in taxids.values()
    ]

    with ThreadPoolExecutor() as executor:
        future_to_task = {
            executor.submit(genbank.fetch_gb_records, *task, count=True): task
            for task in tasks
        }

    results = {}
    errors = []
    for future in as_completed(future_to_task):
        taxid, _ = future_to_task[future]
        try:
            results[taxid] = future.result()
        except Exception as exc:
            logger.error(
                f"[{MODULE_NAME}]: Error processing fetch_gb_records for"
                f" taxid {taxid}:\n{exc}")
            errors.append((taxid_to_species[taxid], exc))

    species_counts = {
        taxid_to_species[taxid]: count
        for taxid, count in results.items()
    }
    # TODO: potential for caching related species counts here
    return species_counts, errors


def write_db_coverage(query_dir, results):
    pass


if __name__ == '__main__':
    main()
