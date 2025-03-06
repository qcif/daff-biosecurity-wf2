"""Assess database coverage of target species at the given locus.

This module features a lot of threading for requests against the Entrez and
GBIF APIs. There can be thousands of requests made in a single run, so a lock
file is used to limit the number of concurrent requests.
"""

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from pprint import pformat

from src.entrez import genbank
from src.utils.flags import Flag, FLAGS, TARGETS
from src.gbif.relatives import GBIFRecordNotFound, RANK, RelatedTaxaGBIF
from src.taxonomy import extract
from src.utils import errors
from src.utils.config import Config

logger = logging.getLogger(__name__)
config = Config()

MODULE_NAME = "Database Coverage"


def read_candidate_species(query_dir):
    candidates = config.read_json(query_dir / config.CANDIDATES_JSON)
    return [
        c["species"]
        for c in candidates["species"]
    ]


def _get_targets(query_dir):
    candidates = read_candidate_species(query_dir)
    if len(candidates) > config.DB_COVERAGE_MAX_CANDIDATES:
        logger.info(
            f"[{MODULE_NAME}]: Skipping database coverage assessment for"
            f" candidates: more than {config.DB_COVERAGE_MAX_CANDIDATES}"
            f" candidates species have been identified ({len(candidates)})."
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
    return candidates, toi_list, pmi


def _get_taxids(targets, query_dir):
    target_taxids = extract.taxids(targets)
    if len(target_taxids) != len(targets):
        msg = (
            "Taxonkit failed to produce taxids for some target species."
            " Database coverage for these species is assumed to be zero, since"
            " this likely means they are not represented in the reference"
            " database.")
        logger.warning(
            f"[{MODULE_NAME}]: {msg}")
        errors.write(
            errors.LOCATIONS.DATABASE_COVERAGE_NO_TAXID,
            msg,
            None,
            query_dir=query_dir,
            data={
                "targets": [k for k, v in target_taxids.items() if not v],
            },
        )
    return target_taxids


def _fetch_target_taxa(targets, query_dir):
    target_gbif_taxa = {}
    higher_taxon_targets = {}  # Taxa at rank 'family' or higher
    for target in targets:
        try:
            gbif_target = RelatedTaxaGBIF(target)
        except GBIFRecordNotFound as exc:
            msg = (f"No GBIF record found for target species '{target}'."
                   " This target could not be evaluated.")
            logger.warning(
                f"[{MODULE_NAME}]: {msg}")
            errors.write(
                errors.LOCATIONS.DATABASE_COVERAGE_NO_GBIF_RECORD,
                msg,
                exc,
                query_dir=query_dir,
                data={"target": target},
            )
            continue
        if gbif_target.rank > RANK.GENUS:
            # These get processed differently - broad GB record count only
            higher_taxon_targets[target] = gbif_target
        else:
            target_gbif_taxa[target] = gbif_target

    logger.debug(
        f"[{MODULE_NAME}]: Targets identified at rank genus or lower:\n"
        + pformat(list(target_gbif_taxa.keys()), indent=2)
    )
    logger.debug(
        f"[{MODULE_NAME}]: Targets identified at rank family or higher:\n"
        + pformat(list(higher_taxon_targets.keys()), indent=2)
    )

    return target_gbif_taxa, higher_taxon_targets


def _parallel_process_tasks(
    tasks,
    query_dir,
    target_taxids,
    target_gbif_taxa,
    taxid_to_taxon,
    candidate_list,
    toi_list,
    pmi,
):
    with ThreadPoolExecutor(max_workers=50) as executor:
        results = {
            get_target_coverage.__name__: {},
            get_related_coverage.__name__: {},
            get_related_country_coverage.__name__: {},
        }
        logger.debug(
            f"[{MODULE_NAME}]: Threading {len(tasks)} tasks..."
        )
        future_to_task = {
            executor.submit(*task): task
            for task in tasks
        }
        for future in as_completed(future_to_task):
            func, target = future_to_task[future][:2]
            logger.info(f"Task completed: {func.__name__} on target"
                        f" '{target}'")
            try:
                results[func.__name__][target] = future.result()
            except Exception as exc:
                species_name = (
                    taxid_to_taxon[target]
                    if isinstance(target, str)
                    else target
                )
                target_source = (
                    "candidate" if species_name in candidate_list
                    else "taxon of interest" if species_name in toi_list
                    else "preliminary ID"
                )
                msg = (
                    f"Error processing {func.__name__} for target species"
                    f" '{species_name}' ({target_source}). This target could"
                    f" not be evaluated: {exc}")
                logger.error(f"[{MODULE_NAME}]: {msg}")
                errors.write(
                    errors.LOCATIONS.DATABASE_COVERAGE,
                    msg,
                    exc,
                    query_dir=query_dir)
                if 'failure in name resolution' in str(exc):
                    raise errors.APIError(
                        f"Fatal error fetching data from Entrez API: '{exc}'"
                        " This error occurred multiple times and indicates a"
                        " network issue - please resume this"
                        " job at a later time when network issues have"
                        " resolved. If this issue persists, you may need to"
                        " contact the development team to diagnose the"
                        " issue. You can check the status of the Entrez API"
                        " by visiting"
                        " https://eutils.ncbi.nlm.nih.gov"
                        "/entrez/eutils/efetch.fcgi in your browser.")

    logger.debug("Results collected from tasks:")
    for func, result in results.items():
        for k in result:
            logger.debug(f"{func}: {k}")

    return _collect_results(
        results,
        target_taxids,
        target_gbif_taxa,
        taxid_to_taxon,
        candidate_list,
        toi_list,
        pmi,
    )


def _collect_results(
    results,
    target_taxids,
    target_gbif_taxa,
    taxid_to_taxon,
    candidate_list,
    toi_list,
    pmi,
):
    error_detected = False
    candidate_results = {}
    toi_results = {}
    pmi_results = {}

    for taxid in target_taxids.values():
        target_taxon = taxid_to_taxon[taxid]
        result = results[get_target_coverage.__name__].get(taxid)
        error_detected = error_detected or result is None

        if target_taxon in candidate_list:
            candidate_results[target_taxon] = candidate_results.get(
                target_taxon, {})
            candidate_results[target_taxon]['target'] = result
        if target_taxon in toi_list:
            toi_results[target_taxon] = toi_results.get(target_taxon, {})
            toi_results[target_taxon]['target'] = result
        if target_taxon == pmi:
            pmi_results[target_taxon] = {
                'target': result,
            }

    for target_taxon, gbif_taxon in target_gbif_taxa.items():
        related_result = results[get_related_coverage.__name__].get(gbif_taxon)
        country_result = results[get_related_country_coverage.__name__].get(
            gbif_taxon)
        error_detected = (
            error_detected
            or related_result is None
            or country_result is None
        )
        if target_taxon in candidate_list:
            candidate_results[target_taxon]['related'] = related_result
            candidate_results[target_taxon]['country'] = country_result
        if target_taxon in toi_list:
            toi_results[target_taxon]['related'] = related_result
            toi_results[target_taxon]['country'] = country_result
        if target_taxon == pmi:
            pmi_results[target_taxon]['related'] = related_result
            pmi_results[target_taxon]['country'] = country_result

    return {
        TARGETS.CANDIDATE: candidate_results,
        TARGETS.TOI: toi_results,
        TARGETS.PMI: pmi_results,
    }, error_detected


def assess_coverage(query_dir) -> dict[str, dict[str, dict]]:
    def get_args(func, query_dir, target, taxid, locus, country):
        if func == get_target_coverage:
            return func, taxid, locus
        elif func == get_related_coverage:
            return func, target, locus, query_dir
        elif func == get_related_country_coverage:
            return func, target, locus, country, query_dir

    locus = config.get_locus_for_query(query_dir)
    country = config.get_country_code_for_query(query_dir)
    candidate_list, toi_list, pmi = _get_targets(query_dir)
    targets = candidate_list + toi_list + [pmi]
    if not targets:
        logger.info(
            f"[{MODULE_NAME}]: Skipping analysis - no target species"
            " identified for database coverage assessment."
        )
        return None

    target_taxids = _get_taxids(targets, query_dir)
    taxid_to_taxon = {v: k for k, v in target_taxids.items()}
    logger.info(
        f"[{MODULE_NAME}]: Assessing database coverage for {len(targets)}"
        f" species at locus '{locus}' in country '{country}'."
    )
    logger.debug(
        f"[{MODULE_NAME}]: collected targets:\n"
        f"  - Candidates: {candidate_list}\n"
        f"  - Taxa of interest: {toi_list}\n"
        f"  - PMI: {pmi}"
    )
    logger.debug(
        f"[{MODULE_NAME}]: Taxids for targets (extracted by taxonkit):\n"
        + pformat(target_taxids, indent=2)
    )

    # 'Higher taxa' are at rank 'family' or higher
    target_gbif_taxa, higher_taxon_targets = _fetch_target_taxa(
        targets, query_dir)

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
        if target in target_gbif_taxa
    ]

    tasks += [
        (get_target_coverage, taxid, locus)
        for target, taxid in target_taxids.items()
        if target in higher_taxon_targets
    ]

    if not len(tasks):
        raise ValueError(
            "No tasks created for database coverage assessment. This likely"
            " indicates a bug in the code - please report this issue.")

    results, is_error = _parallel_process_tasks(
        tasks,
        query_dir,
        target_taxids,
        target_gbif_taxa,
        taxid_to_taxon,
        candidate_list,
        toi_list,
        pmi,
    )
    _set_flags(results, query_dir)
    return results, is_error


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
        for r in gbif_target.relatives
    })
    if not species_names:
        return {}
    logger.info(
        f"[{MODULE_NAME}]: Fetching Genbank records for target"
        f" '{gbif_target.taxon}' (locus: '{locus}') - {len(species_names)}"
        f" related species..."
    )
    # TODO: potential for caching GBIF related taxa here
    results, err = fetch_gb_records_for_species(species_names, locus)
    if err:
        for species, exc in err:
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
        return {}
    # TODO: potential for caching GBIF related/country taxa here
    logger.info(
        f"[{MODULE_NAME}]: Fetching Genbank records for target"
        f" '{gbif_target.taxon}' (locus: '{locus}'; country: '{country}')"
        f" - {len(species_names)} related species"
    )
    results, err = fetch_gb_records_for_species(species_names, locus)
    if err:
        for species, exc in err:
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


def _set_flags(db_coverage, query_dir):
    """Set flags 5.1 - 5.3 (DB coverage) for each target."""
    def set_target_coverage_flag(target, target_type, count):
        if count is None:
            # TODO: this would indicate an error which should be reported
            # elsewhere
            return
        if count > config.CRITERIA.DB_COV_TARGET_MIN_A:
            flag_value = FLAGS.A
        elif count > config.CRITERIA.DB_COV_TARGET_MIN_B:
            flag_value = FLAGS.B
        else:
            flag_value = FLAGS.C
        Flag.write(
            query_dir,
            FLAGS.DB_COVERAGE_TARGET,
            flag_value,
            target=target,
            target_type=target_type,
        )

    def set_related_coverage_flag(target, target_type, species_counts):
        if species_counts is None:
            return  # TODO: Indicates a fatal error
        if not species_counts:
            return  # TODO: no species to check? Flag D??
        total_species = len(species_counts)
        represented_species = len([
            count for count in species_counts.values()
            if count > 0
        ])
        percent_coverage = 100 * represented_species / total_species
        if percent_coverage > config.CRITERIA.DB_COV_RELATED_MIN_A:
            flag_value = FLAGS.A
        elif percent_coverage > config.CRITERIA.DB_COV_RELATED_MIN_B:
            flag_value = FLAGS.B
        else:
            flag_value = FLAGS.C
        Flag.write(
            query_dir,
            FLAGS.DB_COVERAGE_RELATED,
            flag_value,
            target=target,
            target_type=target_type,
        )

    def set_country_coverage_flag(target, target_type, species_counts):
        if species_counts is None:
            return  # TODO: Indicates a fatal error
        total_species = len(species_counts)
        represented_species = len([
            count for count in species_counts.values()
            if count > 0
        ])
        unrepresented_species = total_species - represented_species
        if not species_counts:
            flag_value = FLAGS.C
        elif unrepresented_species <= config.CRITERIA.DB_COV_COUNTRY_MISSING_A:
            flag_value = FLAGS.A
        elif unrepresented_species <= config.CRITERIA.DB_COV_COUNTRY_MISSING_B:
            flag_value = FLAGS.B
        Flag.write(
            query_dir,
            FLAGS.DB_COVERAGE_RELATED_COUNTRY,
            flag_value,
            target=target,
            target_type=target_type,
        )

    for target_type, target_data in db_coverage.items():
        for target_species, coverage_data in target_data.items():
            if not coverage_data.get('related'):
                # No related coverage - target rank is family or higher
                continue
            set_target_coverage_flag(
                target_species,
                target_type,
                coverage_data['target'],
            )
            set_related_coverage_flag(
                target_species,
                target_type,
                coverage_data['related'],
            )
            set_country_coverage_flag(
                target_species,
                target_type,
                coverage_data['country'],
            )
