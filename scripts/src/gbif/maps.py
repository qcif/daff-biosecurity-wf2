from pathlib import Path
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from pygbif import occurrences
import geopandas as gpd
import pandas as pd
from mpl_toolkits.basemap import Basemap
from shapely.geometry import Point
import fsspec

# TAXON_KEY = 2440326
# NAUTRALEARTH_LOWRES_URL = "https://www.naturalearthdata.com/download/110m/cultural/ne_110m_admin_0_countries.zip"
# NAUTRALEARTH_LOWRES_URL = "naturalearth_lowres/ne_110m_admin_0_countries.shp"
NAUTRALEARTH_LOWRES_URL = "scripts/src/gbif/ne_110m_admin_0_countries.zip"


def fetch_gbif_map(taxon_key: str, path: Path):

    res = occurrences.search(taxonKey=taxon_key, limit=1000)
    # results = res["result"]
    lats = [record['decimalLatitude']
            for record in res['results']
            if 'decimalLatitude' in record]
    lons = [record['decimalLongitude']
            for record in res['results']
            if 'decimalLongitude' in record]

    df = pd.DataFrame({'latitude': lats, 'longitude': lons})
    # geometry = [Point(xy) for xy in zip(df['longitude'], df['latitude'])]
    # gdf_points = gpd.GeoDataFrame(df, geometry=geometry, crs="EPSG:4326")

    # df = pd.DataFrame(res['results'])

    # # Convert to GeoDataFrame
    gdf_points = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['longitude'], df['latitude']))

    # # Plot occurrences
    with fsspec.open(f"simplecache::{NAUTRALEARTH_LOWRES_URL}") as file:
        world = gpd.read_file(file)
    fig, ax = plt.subplots(figsize=(14, 10))
    world.plot(ax=ax, color='lightgray')  # Plot the world map
    # Create a hexbin plot to mimic GBIF density visualization
    hb = ax.hexbin(
        x=df['longitude'], y=df['latitude'], gridsize=50,
        cmap='Reds', mincnt=1, alpha=0.8, norm=LogNorm()  # Adjust color and transparency
    )

    # Add a colorbar to show density scale
    cb = fig.colorbar(hb, ax=ax, orientation='vertical')
    cb.set_label("Density of Occurrences")

    # Add title and show the map
    plt.title("Density Map of Species Occurrences (Taxon Key: {taxonKey})")
    plt.show()


if __name__ == "__main__":
    # Example usage
    TAXON_KEY = "2440326"
    OUTPUT_PATH = Path("output/density_map_species.png")

    # Call the function to fetch and save the map
    fetch_gbif_map(taxon_key=TAXON_KEY, path=OUTPUT_PATH)
