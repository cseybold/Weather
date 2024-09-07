'''import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Open the grib files
ds1 = xr.open_dataset('panguweatherCreek.grib', engine='cfgrib')
ds2 = xr.open_dataset('fourcastnetCreek.grib', engine='cfgrib')

# Select the data at 1000 hPa pressure level
temp1 = ds1['t'].sel(isobaricInhPa=850)
temp2 = ds2['t'].sel(isobaricInhPa=850)

# Interpolate temp2 to the coordinates of temp1
temp2 = temp2.interp_like(temp1)

# Calculate the difference
diff = temp1 - temp2

# Select the first time step
diff = diff.isel(step=11)

print(diff.min().values, diff.max().values)
print(ds1.valid_time.values)
print(ds2.keys())

# Create a figure
fig = plt.figure(figsize=(10, 8))

# Create a map projection
ax = plt.axes(projection=ccrs.PlateCarree())

# Add coastlines
ax.coastlines()

# Add country boundaries
ax.add_feature(cfeature.BORDERS, edgecolor='black', alpha=0.5)

# Add state boundaries
states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none')
ax.add_feature(states_provinces, edgecolor='black', alpha=0.5)

# Plot the data
img = diff.plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), cmap='coolwarm', add_labels=False)

# Add a colorbar
plt.colorbar(img, ax=ax, orientation='vertical', label='Temperature Difference (K)')

# Add a title to the plot
plt.title('Temperature Difference (K) at 1000 hPa', fontsize=16)

# Show the plot
plt.show()
'''

import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.gridspec import GridSpec

# Open the grib files
ds1 = xr.open_dataset('panguweatherCreek.grib', engine='cfgrib')
ds2 = xr.open_dataset('fourcastnetCreek.grib', engine='cfgrib')

# Select the data at 850 hPa pressure level
temp1 = ds1['t'].sel(isobaricInhPa=850)
temp2 = ds2['t'].sel(isobaricInhPa=850)

# Interpolate temp2 to the coordinates of temp1
temp2 = temp2.interp_like(temp1)

# Calculate the difference
diff = temp1 - temp2

# Create a figure with 3 subplots using GridSpec
fig = plt.figure(figsize=(12, 18))
gs = GridSpec(3, 2, width_ratios=[1, 0.05], height_ratios=[1, 1, 1], wspace=0.1, hspace=0.3)

# Define the steps to plot
steps = [3, 7, 11]

# Titles for the subplots
titles = ['Temperature Difference (K) After Day 1', 'Temperature Difference (K) After Day 2',
          'Temperature Difference (K) After Day 3']

for i, step in enumerate(steps):
    ax = fig.add_subplot(gs[i, 0], projection=ccrs.PlateCarree())
    diff_step = diff.isel(step=step)

    # Add coastlines
    ax.coastlines()

    # Add country boundaries
    ax.add_feature(cfeature.BORDERS, edgecolor='black', alpha=0.5)

    # Add state boundaries
    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')
    ax.add_feature(states_provinces, edgecolor='black', alpha=0.5)

    # Plot the data
    img = diff_step.plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), cmap='coolwarm', add_labels=False,
                                    add_colorbar=False)

    # Add a title to the subplot
    ax.set_title(titles[i], fontsize=16)

# Add a single colorbar for the entire figure
cbar_ax = fig.add_subplot(gs[:, 1])
cbar = fig.colorbar(img, cax=cbar_ax, orientation='vertical', label='Temperature Difference (K)')

# Adjust layout
plt.subplots_adjust(left=0.05, right=0.9, top=0.95, bottom=0.05)

# Show the plot
plt.show()

