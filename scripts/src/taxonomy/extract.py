import logging
import subprocess
import tempfile
from pathlib import Path

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


def taxonomies(
    taxids: list[str],
    taxdb: Path = Path(config.TAXONKIT_DATA),
) -> dict[str, dict[str, str]]:
    """Use taxonkit lineage to extract taxonomic data for given taxids."""
    with tempfile.NamedTemporaryFile(mode='w+') as temp_file:
        temp_file.write("\n".join(taxids))
        temp_file.flush()
        temp_file_name = temp_file.name
        result = subprocess.run(
            [
                'taxonkit',
                'lineage',
                '-R',
                '-c', temp_file_name,
                '--data-dir', taxdb,
            ],
            capture_output=True,
            text=True
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
    with tempfile.NamedTemporaryFile(mode='w+') as temp_file:
        temp_file.write("\n".join(species_list))
        temp_file.flush()
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

    logger.debug(
        "[extract.taxids] taxonkit name2taxid stdout returned:\n"
        + result.stdout.strip()
    )
    if result.stderr.strip():
        logger.warning(
            "[extract.taxids] taxonkit name2taxid stderr returned:\n"
            + result.stderr.strip()
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
