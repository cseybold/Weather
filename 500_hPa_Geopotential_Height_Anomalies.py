import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import math
import cartopy.feature as cfeature

# Open the grib file
ds = xr.open_dataset('panguweatherPaper72.grib', engine='cfgrib')

# Select the variable you want to plot (change 'z' to your variable)
da = ds['z']

# Select a specific pressure level (e.g., 500 hPa)
da = da.sel(isobaricInhPa=500)

# Compute anomalies in standard deviations
mean = da.mean(dim='step')
std_dev = da.std(dim='step')
da = (da - mean) / std_dev

# Get the number of time steps
num_steps = len(da.step) - 1

# Number of days and timesteps per day
timesteps_per_day = 4
num_days = math.ceil(num_steps/timesteps_per_day)

fig, axs = plt.subplots(num_days, timesteps_per_day, figsize=(15, 15), subplot_kw={'projection': ccrs.PlateCarree()})

# Define the latitude and longitude boundaries for Western North America
lon_min, lon_max = -150, -100
lat_min, lat_max = 20, 70

# Create a subplot for each time step
for i in range(num_steps):
    # Calculate the current row and column
    col = i % timesteps_per_day
    row = i // timesteps_per_day

    # Select the data at this time step
    da_step = da.isel(step=i)

    # Add coastlines to the map
    axs[row, col].coastlines()

    # Set the extent to the boundaries of Western North America
    axs[row, col].set_extent([lon_min, lon_max, lat_min, lat_max])

    # Add country boundaries
    axs[row, col].add_feature(cfeature.BORDERS, edgecolor='black', alpha=0.5)

    # Add state boundaries
    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')
    axs[row, col].add_feature(states_provinces, edgecolor='black', alpha=0.5)

    # Add gridlines
    gl = axs[row, col].gridlines(draw_labels=True, linewidth=1, color='gray', alpha=0.5)
                                 #, xlocs=np.arange(lon_min, lon_max, 10), ylocs=np.arange(lat_min, lat_max, 10))
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = False
    gl.left_labels = False

    # Plot the data
    img = da_step.plot(ax=axs[row, col], transform=ccrs.PlateCarree(), add_colorbar=False, add_labels=False)

    # Add time labels to the rows and day labels to the columns
    if row == 0:
        time_str = str((col % timesteps_per_day) * 6).zfill(2) + " UTC"
        axs[row, col].set_title(time_str, fontsize=20, fontweight='bold')

    if col == 0:
        date_str = da_step.valid_time.values.astype('datetime64[D]').astype(str)
        fig.text(0.1, 0.72 - (row / timesteps_per_day) * 1.01, date_str,
                 rotation='vertical', fontsize=20, fontweight='bold')

    if row == 2:
        gl.bottom_labels = True
    if col == 3:
        gl.right_labels = True

# Remove unused subplots
for i in range(num_steps, num_days*timesteps_per_day):
    fig.delaxes(axs.flatten()[i])

# Add a single colorbar to the figure
cbar_ax = fig.add_axes([0.1, 0.05, 0.8, 0.02])  # This creates a new axes where the colorbar will be
cbar = fig.colorbar(img, cax=cbar_ax, orientation='horizontal')  # This adds the colorbar to the newly created axes
cbar.set_label('σ', size=15, fontweight='bold')  # Replace 'Your Label' with your actual label

# Adjust the layout to prevent overlap
fig.subplots_adjust(bottom=0.15)

# Add a title to the plot
plt.suptitle('500 hPa Geopotential Height Anomalies', fontsize=32, fontweight='bold')

# Show the plot
plt.show()
