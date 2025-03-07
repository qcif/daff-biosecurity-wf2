from pathlib import Path
from matplotlib import pyplot as plt
from pygbif import occurrences
import geopandas as gpd
import pandas as pd
import fsspec

NATURALEARTH_LOWRES_URL = "scripts/src/gbif/ne_110m_admin_0_countries.zip"


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
    gdf = gpd.GeoDataFrame(df,
                           geometry=gpd.points_from_xy(df['longitude'],
                                                       df['latitude']))

    with fsspec.open(f"simplecache::{NATURALEARTH_LOWRES_URL}") as file:
        world = gpd.read_file(file)

    ax = world.plot(figsize=(10, 6), color='white', edgecolor='black')
    gdf.plot(ax=ax, color='red', markersize=10, alpha=0.5, label='Occurrences')
    plt.title(f"Species Occurrence in the World (Taxon Key: {taxon_key})")
    plt.legend()
    plt.savefig(path, bbox_inches='tight', dpi=300)
    plt.close()  # Close the plot to free memory


# Example usage
if __name__ == "__main__":
    TAXON_KEY = "2440326"
    OUTPUT_PATH = Path("scripts/src/gbif/species_map.png")
    fetch_gbif_map(taxon_key=TAXON_KEY, path=OUTPUT_PATH)
