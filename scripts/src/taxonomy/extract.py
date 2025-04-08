import logging
import os
import subprocess
import tempfile

from ..utils.config import Config

logger = logging.getLogger(__name__)
config = Config()


TAXONOMIC_RANKS = [
    "superkingdom",
    'kingdom',
    'phylum',
    'class',
    'order',
    'family',
    'genus',
    'species',
]


def taxonomies(taxids: list[str]) -> dict[str, dict[str, str]]:
    """Use taxonkit lineage to extract taxonomic data for given taxids."""

    # Because temporary file handling in Windows is different,
    # delete parameter need to be set to False and closed manually
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_file:
        temp_file.write("\n".join(taxids))
        temp_file.flush()
        temp_file_name = temp_file.name
    try:
        result = subprocess.run(
            [
                'taxonkit',
                'lineage',
                '-R',
                '-c', temp_file_name,
                '--data-dir', config.TAXONKIT_DATA,
            ],
            capture_output=True,
            check=True,
            text=True,
        )
    except subprocess.CalledProcessError as exc:
        logger.error(
            "[extract.taxonomies] taxonkit lineage failed with error:\n"
            + exc.stderr
        )
        raise exc
    finally:
        if temp_file:
            temp_file.close()
        # Clean up the temporary file after processing
        if os.path.exists(temp_file_name):
            os.remove(temp_file_name)
            logger.info(
                f"Temporary file {temp_file_name} deleted successfully."
            )

    logger.debug(
        "[extract.taxids] taxonkit name2taxid stdout:\n"
        + result.stdout
    )
    if result.stderr.strip():
        logger.warning(
            "[extract.taxids] taxonkit name2taxid stderr:\n"
            + result.stderr
        )

    taxonomy_data = {}
    for line in result.stdout.strip().split('\n'):
        fields = line.split('\t')[1:]
        if len(fields) == 3:
            taxid, taxon_details, ranks = fields[0], fields[1], fields[2]
            lineage_list = taxon_details.split(';')
            ranks_list = ranks.split(';')
            taxonomy = {
                rank: name for rank,
                name in zip(ranks_list, lineage_list)
                if rank in TAXONOMIC_RANKS
            }
            taxonomy_data[taxid] = taxonomy
        else:
            logger.warning(
                "[extract.taxonomies] Warning: Unexpected format in taxonkit"
                " stdout. This may result in missing taxonomy information:\n"
                + line)
    return taxonomy_data


def taxids(species_list: list[str]) -> dict[str, str]:
    """Use taxonkit name2taxid to extract taxids for given species.

    These species did not come from the core_nt database, so they might not
    even have a taxid if they are unsequenced/rare/new species.
    """
    logger.debug(
        "[extract.taxids] Extracting taxids for species using taxonkit"
        f" name2taxid with {len(species_list)} species:\n"
        + "\n".join(species_list[:3] + ['...'])
    )

    # Because temporary file handling in Windows is different,
    # delete parameter need to be set to False and closed manually
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_file:
        temp_file.write("\n".join(species_list))
        temp_file.flush()
        temp_file_name = temp_file.name
    try:
        result = subprocess.run(
            [
                'taxonkit',
                'name2taxid',
                temp_file.name,
                '--data-dir', config.TAXONKIT_DATA,
            ],
            capture_output=True,
            text=True,
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        logger.error(
            "[extract.taxids] taxonkit name2taxid failed with error:\n"
            + exc.stderr
        )
        raise exc
    finally:
        if temp_file:
            temp_file.close()
        # Clean up the temporary file after processing
        if os.path.exists(temp_file_name):
            os.remove(temp_file_name)
            logger.info(
                f"Temporary file {temp_file_name} deleted successfully."
            )

    logger.debug(
        "[extract.taxids] taxonkit name2taxid stdout:\n"
        + result.stdout
    )
    if result.stderr.strip():
        logger.warning(
            "[extract.taxids] taxonkit name2taxid stderr:\n"
            + result.stderr
        )

    taxid_data = {}
    for line in result.stdout.strip().split('\n'):
        if not line.strip():
            continue
        fields = line.split('\t')
        if len(fields) == 2:
            species, taxid = fields
            taxid_data[species] = taxid or None
        elif (field := fields[0].strip()) and len(fields) == 1:
            species = field
            taxid_data[species] = None
        else:
            logger.warning(
                "[extract.taxids] Warning: Unexpected format in taxonkit"
                " stdout. This may result in missing taxid information:\n"
                + line)
    return taxid_data
