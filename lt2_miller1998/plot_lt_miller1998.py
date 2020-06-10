
# coding: utf-8

# ## To do
# 
# - Agree on a title naming approach (e.g. do we keep "Laugh test X:", do we want to include the main reference in the title)
# - Check the expectations for this laugh test and update the image accordingly
# - Agree on what constitutes relevant meta data for these things
# - Update results with Reza's new ones
# - Do we want the code to be more elegant in terms of extracting length of data series, convert time from the .nc attributes etc?
# - Image format and specs?

# # Laugh test 2: Storage and transmission in soils
# This notebook plots SUMMA simulations for the test case defined in Miller et al. (1998). 
# 
# ## Expectations
# The simulations should reproduce Figure 1 in Miller et al. (1998):
# 
# <div>
# <img src="img/miller1998_fig1.png" width="300"/>
#     <center> Figure 1a in Miller et al. (1998): "Dense-grid solutions" </center>
# </div>
# 
# ## Workflow
# - Load model simulations (.nc) into memory
# - Extract `mLayerMatricHead` and `mLayerHeight` variables for experiments with sand, loam and clay
# - Create a plot of pressure head vs depth, for the selected times
# 
# ## Meta data
# 
# | Data  | Value  |
# |:---|:---|
# | Model name| Structure for Unifying Multiple Modelling Alternatives (SUMMA) |
# | Model version  |   |
# | Model reference | Clark et al. (2015a,b) |
# | Model runs by | M. Clark |
# | Model decisions | see attributes in output.nc  |
# | Model run date | see attributes in output.nc |
# | Notebook code by | W. Knoben, A. Bennett |
#     
# ## Reference(s)
# Miller, C. T., Williams, G. A., Kelley, C. T. & Tocci, M. D. (1998), Robust solution of Richard's equation for nonuniform porous media, Water Resour. Res., 34(10), 2599–2610, https://doi.org/10.1029/98WR01673
# 
# Clark, M. P., Nijssen, B., Lundquist, J. D., Kavetski, D., Rupp, D. E., Woods, R. A., … Rasmussen, R. M. (2015a). A unified approach for process-based hydrologic modeling: 1. Modeling concept. Water Resources Research, 51(4), 2498–2514. https://doi.org/10.1002/2015WR017198
# 
# Clark, M. P., Nijssen, B., Lundquist, J. D., Kavetski, D., Rupp, D. E., Woods, R. A., … Marks, D. G. (2015b). A unified approach for process-based hydrologic modeling: 2. Model implementation and case studies. Water Resources Research, 51, 2515–2542. https://doi.org/10.1002/2015WR017200

# In[80]:


# modules
from pathlib import Path
from matplotlib.ticker import MultipleLocator
import xarray as xr # note, also needs netcdf4 library installed
import matplotlib.pyplot as plt


# In[2]:


# Specify the data locations relative to the notebook
sim_path = Path("./output")
sim_sand = "millerSand_output_timestep.nc"
sim_clay = "millerClay_output_timestep.nc"
sim_loam = "millerLoam_output_timestep.nc"


# In[3]:


# Load the data
ds_sand = xr.open_dataset( sim_path / sim_sand )
ds_clay = xr.open_dataset( sim_path / sim_clay )
ds_loam = xr.open_dataset( sim_path / sim_loam )


# In[4]:


# Extract the variables we want
matricHead_sand = ds_sand.mLayerMatricHead
matricHead_clay = ds_clay.mLayerMatricHead
matricHead_loam = ds_loam.mLayerMatricHead

layerHeight_sand = ds_sand.mLayerHeight
layerHeight_clay = ds_clay.mLayerHeight
layerHeight_loam = ds_loam.mLayerHeight


# In[27]:


# Define the time of interest (see test case description), converted into a timestep index by dividing by time step length
dt = 900 #[s]
time_sand = int( (0.18*86400)/dt ) # d * (s/d) / (s/timestep)
time_clay = int( (1.00*86400)/dt ) # int() rounds down to nearest integer
time_loam = int( (2.25*86400)/dt )


# In[45]:


# Select the time slices of interest
matricHead_sand_t17 = matricHead_sand.isel(time=time_sand)
matricHead_clay_t17 = matricHead_clay.isel(time=time_clay)
matricHead_loam_t17 = matricHead_loam.isel(time=time_loam)

layerHeight_sand_t17 = layerHeight_sand.isel(time=time_sand)
layerHeight_clay_t17 = layerHeight_clay.isel(time=time_clay)
layerHeight_loam_t17 = layerHeight_loam.isel(time=time_loam)


# In[46]:


# Reverse SUMMA's layerHeight_ variables 
# These currently show a depth below surface, with positive values for deeper layers)
# We want them to be a height above bottom of the profile (see test case specifications)
layerHeight_sand_t17 = layerHeight_sand_t17[::-1]
layerHeight_clay_t17 = layerHeight_clay_t17[::-1]
layerHeight_loam_t17 = layerHeight_loam_t17[::-1]


# In[29]:


# Ensure that we can actually read the figure labels
font = {'weight' : 'normal',
        'size'   : 18}

plt.rc('font', **font)


# In[81]:


# Open a figure
fig = plt.figure(figsize=(18, 10), dpi= 80, facecolor='w', edgecolor='k');

# Plot the data in a single figure
plt.plot(matricHead_sand_t17[0:800], layerHeight_sand_t17[5:805], marker='', color='blue',   label='sand')
plt.plot(matricHead_loam_t17[0:400], layerHeight_loam_t17[5:405], marker='', color='orange', label='loam')
plt.plot(matricHead_clay_t17[0:320], layerHeight_clay_t17[5:325], marker='', color='green',  label='clay loam')

# Labels
plt.xlabel('Pressure head [m]'); # note, ';' supresses output from the Text object that is created for the labels
plt.ylabel('Height above profile bottom[m]');
plt.legend();

# Axes
plt.xlim(left=-5)
plt.ylim(0,10)
ax = plt.gca()
ax.yaxis.set_minor_locator(MultipleLocator(0.5))
ax.tick_params(which='major', length=7, direction='in')
ax.tick_params(which='minor', length=4, direction='in')

# Save the figure
plt.savefig('img/lt2_miller1998.png');

