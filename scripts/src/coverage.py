"""Assess database coverage of target species at the given locus.

This module features a lot of threading for requests against the Entrez and
GBIF APIs. There can be thousands of requests made in a single run, so a lock
file is used to limit the number of concurrent requests.
"""

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from pprint import pformat

from src.entrez import genbank
from src.utils import errors
from src.utils.config import Config
from src.gbif.relatives import NoGenusFoundError, RelatedTaxaGBIF
from src.taxonomy import extract

logger = logging.getLogger(__name__)
config = Config()

MODULE_NAME = "Database Coverage"


def read_candidate_species(query_dir):
    candidates = config.read_json(query_dir / config.CANDIDATES_JSON)
    return [
        c["species"]
        for c in candidates["species"]
    ]


def assess_coverage(query_dir):
    def get_args(func, query_dir, target, taxid, locus, country):
        if func == get_target_coverage:
            return func, taxid, locus
        elif func == get_related_coverage:
            return func, target, locus, query_dir
        elif func == get_related_country_coverage:
            return func, target, locus, country, query_dir

    locus = config.get_locus_for_query(query_dir)
    country = config.get_country_code_for_query(query_dir)
    candidates = read_candidate_species(query_dir)
    if len(candidates) > config.DB_COVERAGE_MAX_CANDIDATES:
        logger.info(
            f"[{MODULE_NAME}]: Skipping database coverage assessment:"
            f" more than {config.DB_COVERAGE_MAX_CANDIDATES} candidates"
            f" species have been identified ({len(candidates)})."
        )
        candidates = []
    pmi = config.get_pmi_for_query(query_dir)
    toi_list = config.get_toi_list_for_query(query_dir)
    if len(toi_list) > config.DB_COVERAGE_TOI_LIMIT:
        toi_list = toi_list[:config.DB_COVERAGE_TOI_LIMIT]
        excluded_tois = toi_list[config.DB_COVERAGE_TOI_LIMIT:]
        msg = (
            f"Only the first {config.DB_COVERAGE_TOI_LIMIT} taxa of interest"
            f" will be evaluated. The following taxa of interest will be"
            f" excluded: {', '.join(excluded_tois)}. This limit can be raised"
            f" by setting the 'DB_COVERAGE_TOI_LIMIT' environment variable.")
        logger.warning(f"[{MODULE_NAME}]: {msg}")
        errors.write(
            errors.LOCATIONS.DATABASE_COVERAGE,
            msg,
            None,
            query_dir=query_dir,
        )
    targets = candidates + toi_list + [pmi]
    target_taxids = extract.taxids(targets)
    taxid_to_species = {v: k for k, v in target_taxids.items()}

    logger.info(
        f"[{MODULE_NAME}]: Assessing database coverage for {len(targets)}"
        f" species at locus '{locus}' in country '{country}'."
    )

    target_gbif_taxa = {}
    for target in targets:
        try:
            gbif_target = RelatedTaxaGBIF(target)
        except NoGenusFoundError as exc:
            msg = (f"No genus found for target species '{target}'. This target"
                   " could not be evaluated.")
            logger.warning(
                f"[{MODULE_NAME}]: {msg}")
            errors.write(
                errors.LOCATIONS.DATABASE_COVERAGE_NO_GENUS,
                msg,
                exc,
                query_dir=query_dir,
                data={"target": target},
            )
            continue
        target_gbif_taxa[target] = gbif_target

    tasks = [
        get_args(
            func,
            query_dir,
            target_gbif_taxa[target],
            taxid,
            locus,
            country,
        )
        for target, taxid in target_taxids.items()
        for func in (
            get_target_coverage,
            get_related_coverage,
            get_related_country_coverage,
        )
    ]

    with ThreadPoolExecutor() as executor:
        logger.debug(
            f"[{MODULE_NAME}]: Threading tasks: {pformat(tasks)}."
        )
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
                if isinstance(target, str)
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
            candidate_results[species_name] = candidate_results.get(
                species_name, {})
            candidate_results[species_name]['target'] = result
        if species_name in toi_list:
            toi_results[species_name] = toi_results.get(species_name, {})
            toi_results[species_name]['target'] = result
        if species_name == pmi:
            pmi_results[species_name] = {
                'target': result,
            }

    for species_name, gbif_taxon in target_gbif_taxa.items():
        related_result = results[get_related_coverage.__name__][gbif_taxon]
        country_result = results[get_related_country_coverage.__name__][
            gbif_taxon]
        if species_name in candidates:
            candidate_results[species_name]['related'] = related_result
            candidate_results[species_name]['country'] = country_result
        if species_name in toi_list:
            toi_results[species_name]['related'] = related_result
            toi_results[species_name]['country'] = country_result
        if species_name == pmi:
            pmi_results[species_name]['related'] = related_result
            pmi_results[species_name]['country'] = country_result

    return {
        "candidates": candidate_results,
        "tois": toi_results,
        "pmi": pmi_results,
    }


def get_target_coverage(taxid, locus):
    """Return a count of the number of accessions for the given target."""
    # TODO: potential for caching gb count result here
    logger.info(
        f"[{MODULE_NAME}]: Fetching Genbank records for target taxid:"
        f"{taxid}, locus: '{locus}'..."
    )
    return genbank.fetch_gb_records(locus, taxid, count=True)


def get_related_coverage(gbif_target, locus, query_dir):
    """Return a count of the number of related species (same genus) and the
    number of species which have at least one accession in the database.
    """
    species_names = list({
        r["canonicalName"]
        for r in gbif_target.related_species
    })
    if not species_names:
        return {}, []
    logger.info(
        f"[{MODULE_NAME}]: Fetching Genbank records for target"
        f" '{gbif_target.taxon}' (locus: '{locus}') - {len(species_names)}"
        f" related species..."
    )
    # TODO: potential for caching GBIF related taxa here
    results, errors = fetch_gb_records_for_species(species_names, locus)
    if errors:
        for species, exc in errors:
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
    if not species_names:
        return {}, []
    # TODO: potential for caching GBIF related/country taxa here
    logger.info(
        f"[{MODULE_NAME}]: Fetching Genbank records for target"
        f" '{gbif_target.taxon}' (locus: '{locus}'; country: '{country}')"
        f" - {len(species_names)} related species"
    )
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
    species_without_taxid = [
        k for k, v in taxids.items()
        if v is None
    ]
    taxid_to_species = {
        v: k
        for k, v in taxids.items()
        if v is not None
    }
    tasks = [
        (locus, taxid)
        for taxid in taxids.values()
        if taxid is not None
    ]

    with ThreadPoolExecutor() as executor:
        future_to_task = {
            executor.submit(genbank.fetch_gb_records, *task, count=True): task
            for task in tasks
        }

    results = {}
    errors = []
    for future in as_completed(future_to_task):
        locus, taxid = future_to_task[future]
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
    species_counts.update({
        species: 0
        for species in species_without_taxid
    })
    # TODO: potential for caching related species counts here
    return species_counts, errors
