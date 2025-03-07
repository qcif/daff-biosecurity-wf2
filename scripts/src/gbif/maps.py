from matplotlib import pyplot as plt
from pygbif import occurrences
import geopandas as gpd
import pandas as pd
from mpl_toolkits.basemap import Basemap
from shapely.geometry import Point
import fsspec

TAXON_KEY = 2440326
NAUTRALEARTH_LOWRES_URL = "scripts/src/gbif/ne_110m_admin_0_countries.zip"



res = occurrences.search(taxonKey=TAXON_KEY, limit=1000)
# results = res["result"]
lats = [record['decimalLatitude']
        for record in res['results']
        if 'decimalLatitude' in record]
lons = [record['decimalLongitude']
        for record in res['results']
        if 'decimalLongitude' in record]

# # Plot the map
# plt.figure(figsize=(14, 10))
# m = Basemap(projection='mill', llcrnrlat=-60,
#             urcrnrlat=90, llcrnrlon=-180,
#             urcrnrlon=180, resolution='c')
# m.drawcoastlines()
# m.drawcountries()
# m.drawmapboundary(fill_color='aqua')
# m.fillcontinents(color='lightgray', lake_color='aqua')
# # x, y = m(longitude, latitude)
# m.scatter(lats, lons, latlon=True,
#           marker='o', color='red',
#           s=10, zorder=5, label="Occurrences")
# # m.scatter(lats, lons, marker="o", color="red", s=10, label="Occurrences")
# plt.title('Species Occurrence Map')
# plt.legend(loc="lower left")
# plt.show()


df = pd.DataFrame({'latitude': lats, 'longitude': lons})
# geometry = [Point(xy) for xy in zip(df['longitude'], df['latitude'])]
# gdf_points = gpd.GeoDataFrame(df, geometry=geometry, crs="EPSG:4326")

# df = pd.DataFrame(res['results'])

# # Convert to GeoDataFrame
gdf_points = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['longitude'], df['latitude']))

# # Plot occurrences
with fsspec.open(f"simplecache::{NAUTRALEARTH_LOWRES_URL}") as file:
    world = gpd.read_file(file)
fig, ax = plt.subplots(figsize=(12, 8))
world.plot(ax=ax, color='lightgray')  # Plot the world map
gdf_points.plot(ax=ax, color='red', markersize=10, label='Occurrences')
# ax = world.plot(figsize=(10, 6), color='white', edgecolor='black')
# gdf.plot(ax=ax, markersize=10, color='red', alpha=0.5)
plt.title("Species Occurrence in the World")
plt.show()

