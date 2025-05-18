"""Fetch GBIF occurrence data and plot on a world map."""

import logging
from pathlib import Path

logging.getLogger('matplotlib').setLevel(logging.WARNING)

import fsspec  # noqa:E402
import geopandas as gpd  # noqa:E402
import pandas as pd  # noqa:E402
from matplotlib import pyplot as plt  # noqa:E402
from matplotlib.colors import LogNorm  # noqa:E402
from pygbif import occurrences  # noqa:E402

from src.utils.config import Config  # noqa:E402

config = Config()
logger = logging.getLogger(__name__)

NATURALEARTH_LOWRES_URL = (
    Path(__file__).parent / 'base_maps/ne_110m_admin_0_countries.zip')


def draw_occurrence_map(taxon_key: str, path: Path):
    '''Fetch GBIF API to get species world map by using taxonomy ID. '''
    all_results = []
    offset = 0
    while True:
        res = occurrences.search(taxonKey=taxon_key, offset=offset)
        all_results.extend(res.get('results', []))
        if res.get('endOfRecords', True):
            break
        offset += res['limit']
        if offset >= config.GBIF_MAX_OCCURRENCE_RECORDS:
            logger.warning(
                "Maximum number of records reached:"
                f" {config.GBIF_MAX_OCCURRENCE_RECORDS}")
            break

    lats = [record['decimalLatitude']
            for record in all_results
            if 'decimalLatitude' in record]
    lons = [record['decimalLongitude']
            for record in all_results
            if 'decimalLongitude' in record]

    if not lats:
        logger.info(f"No occurrence data found for taxon"
                    f" '{taxon_key}' - no map can be drawn.")
        return

    df = pd.DataFrame({'latitude': lats, 'longitude': lons})

    with fsspec.open(f"simplecache::{NATURALEARTH_LOWRES_URL}") as file:
        world = gpd.read_file(file)

    fig, ax = plt.subplots(figsize=(16, 12))
    fig.patch.set_facecolor('#0d1017')
    ax.set_facecolor('#0d1017')
    world.plot(ax=ax, color='#32363c')
    hb = ax.hexbin(
        x=df['longitude'],
        y=df['latitude'],
        gridsize=(150, 30),
        cmap='autumn_r',
        mincnt=1,
        alpha=0.8,
        norm=LogNorm(),
    )

    # Add a colorbar to show density scale
    cb = fig.colorbar(hb, ax=ax, orientation='vertical', shrink=0.7, aspect=60)
    cb.set_label("Density of Occurrences")
    ax.set_axis_off()

    plt.savefig(path, bbox_inches='tight', dpi=300)
    plt.close()
