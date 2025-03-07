from matplotlib import pyplot as plt
from pygbif import occurrences
import geopandas as gpd
import pandas as pd
import requests
from mpl_toolkits.basemap import Basemap

TAXON_KEY = 2440326

res = occurrences.search(taxonKey=TAXON_KEY, limit=1000)
# results = res["result"]
lats = [record['decimalLatitude']
        for record in res['results']
        if 'decimalLatitude' in record]
lons = [record['decimalLongitude']
        for record in res['results']
        if 'decimalLongitude' in record]

# Plot the map
plt.figure(figsize=(14, 10))
m = Basemap(projection='mill', llcrnrlat=-60,
            urcrnrlat=90, llcrnrlon=-180,
            urcrnrlon=180, resolution='c')
m.drawcoastlines()
m.drawcountries()
m.drawmapboundary(fill_color='aqua')
m.fillcontinents(color='lightgray', lake_color='aqua')
# x, y = m(longitude, latitude)
m.scatter(lats, lons, latlon=True,
          marker='o', color='red',
          s=10, zorder=5, label="Occurrences")
# m.scatter(lats, lons, marker="o", color="red", s=10, label="Occurrences")
plt.title('Species Occurrence Map')
plt.legend(loc="lower left")
plt.show()


# # Get occurrence data
# res = occurrences.search(taxonKey=2440326, limit=1000)
# df = pd.DataFrame(res['results'])

# # Convert to GeoDataFrame
# gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['decimalLongitude'], df['decimalLatitude']))

# # Plot occurrences
# world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
# ax = world.plot(figsize=(10, 6), color='white', edgecolor='black')
# gdf.plot(ax=ax, markersize=10, color='red', alpha=0.5)
# plt.title("Species Occurrence in Canada")
# plt.show()

