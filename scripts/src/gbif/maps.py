from pathlib import Path
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from pygbif import occurrences
import geopandas as gpd
import pandas as pd
import fsspec

NATURALEARTH_LOWRES_URL = Path(__file__).parent / 'ne_110m_admin_0_countries.zip'


def fetch_gbif_map(taxon_key: str, path: Path):
    '''Fetch GBIF API to get species world map by using taxonomy ID. '''
    res = occurrences.search(taxonKey=taxon_key, limit=1000)
    lats = [record['decimalLatitude']
            for record in res['results']
            if 'decimalLatitude' in record]
    lons = [record['decimalLongitude']
            for record in res['results']
            if 'decimalLongitude' in record]
    df = pd.DataFrame({'latitude': lats, 'longitude': lons})
    # gdf = gpd.GeoDataFrame(df,
    #                        geometry=gpd.points_from_xy(df['longitude'],
    #                                                    df['latitude']))

    with fsspec.open(f"simplecache::{NATURALEARTH_LOWRES_URL}") as file:
        world = gpd.read_file(file)
    
    fig, ax = plt.subplots(figsize=(16, 12))
    fig.patch.set_facecolor('#d0e8f2')
    ax.set_facecolor('#d0e8f2')
    world.plot(ax=ax, color='lightgray')
    # ax = world.plot(figsize=(10, 6), color='white', edgecolor='black')
    # gdf.plot(ax=ax, color='red', markersize=10, alpha=0.5, label='Occurrences')
    hb = ax.hexbin(
        x=df['longitude'], y=df['latitude'], gridsize=50,
        cmap='Reds', mincnt=1, alpha=0.8, norm=LogNorm()  # Adjust color and transparency
    )

    # Add a colorbar to show density scale
    cb = fig.colorbar(hb, ax=ax, orientation='vertical', shrink=0.7)
    cb.set_label("Density of Occurrences")
    ax.set_axis_off()

    plt.savefig(path, bbox_inches='tight', dpi=300)
    plt.close()  # Close the plot to free memory


# Example usage
if __name__ == "__main__":
    TAXON_KEY = "2440326"
    OUTPUT_PATH = Path("output/query_001_LC438549/species_map.png")
    fetch_gbif_map(taxon_key=TAXON_KEY, path=OUTPUT_PATH)
