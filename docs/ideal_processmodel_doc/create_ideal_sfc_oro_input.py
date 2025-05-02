#!/usr/bin/env python3
# NOTE: requires the netcdf4 package, which is not standard on system installs. 
# On hera, try /home/Ted.Mansell/miniconda3/bin/python3 which has netcdf installed, or 
# otherwise install netcdf in your own miniconda
#
# Purpose: Creates sfc_data.nc and oro_data.nc for UFS ideal setup with user-specified 
# values (affecting land surface, radiation, etc.)
#   THIS WILL NOT OVERWRITE EXISITNG FILES (Error code will be issued)
# nx and ny can be any value greater than or equal to the npy, npy in the fv_core_nml 
#
# The values of lat and lon here don't matter (yet) for ideal cases. Use deglat 
# and deglon in the fv_core_nml namelist to set these (probably only important for
# radiation and Coriolis, if enabled)
#
from netCDF4 import Dataset    # Note: python is case-sensitive!
import numpy as np

# filenames can be anything, but then must link as INPUT/sfc_data.nc and INPUT/oro_data.nc
sfcfilename = 'sfc_data.nc'
orofilename = 'oro_data.nc'

# dimensions:
nx = 300 # should have nx >= npx
ny = 300 # should have ny >= npy
nz = 4  # number of soil layer; should be 4 for noahmp; Can set to 9 for ruc lsm and set namelist appropriately

# The following are used to set various array values 
tmpslmsk = 1 # sea/land mask (0=sea, 1=land)
tmpstype = 1 # land surface type; if set to 0 or 14, will get water surface
# soil types: 1, 'SAND' : 2, 'LOAMY SAND' : 3, 'SANDY LOAM' : 4, 'SILT LOAM' : 5, 'SILT' : 6, 'LOAM' : 7, 'SANDY CLAY LOAM' : 8, 'SILTY CLAY LOAM' : 9, 'CLAY LOAM' : 10,'SANDY CLAY' : 11,'SILTY CLAY' : 12,'CLAY' : 13,'ORGANIC MATERIAL' : 14,'WATER'
tmpvtype = 7 # vegetation type; see noahmptable.tbl; 7 = grassland
tmpvfrac = 0.9 # vegetation fraction
tmpsmc = 0.2 # soil moisture content; set to 1 for water surface; sets all levels to the same value
tmpslc = 0.2 # soil liquid content; set to 0 for water surface; sets all levels to the same value
tmpstc = 300 # soil temperature (all levels); sets all levels to the same value
lat1 = 35  # value is ignored by FV3, which will use deglat from input namelist
lon1 = 270 # value is ignored by FV3, which will use deglon from input namelist
tem2m = 300 # t2m = 2-meter temperature
q2mtmp = 0.009 # 2-meter specific humidity
ust = 0.2 # uustar : for surface layer scheme; cannot be zero

# Note that the variable names (key) are always the same as the long names, so longName is not actually used
sfc_vars = {   'slmsk':  {'longName':'slmsk','units':'none','ndims':'2','value':tmpslmsk}, # sea/land mask array (sea:0,land:1,sea-ice:2)
               'tsea':   {'longName':'tsea','units':'none','ndims':'2','value':tmpstc}, # tsea is read into variable 'tsfc' : surface air temperature
               'sheleg': {'longName':'sheleg','units':'none','ndims':'2','value':0}, # read into variable 'weasd' : water equiv of accumulated snow depth (kg/m**2)
               'tg3':    {'longName':'tg3','units':'none','ndims':'2','value':300}, # deep soil temperature
               'zorl':   {'longName':'zorl','units':'none','ndims':'2','value':0.0001}, # composite surface roughness in cm
               'alvsf':  {'longName':'alvsf','units':'none','ndims':'2','value':0.06}, # mean vis albedo with strong cosz dependency
               'alvwf':  {'longName':'alvwf','units':'none','ndims':'2','value':0.06}, # mean vis albedo with weak cosz dependency
               'alnsf':  {'longName':'alnsf','units':'none','ndims':'2','value':0.06}, # mean nir albedo with strong cosz dependency
               'alnwf':  {'longName':'alnwf','units':'none','ndims':'2','value':0.06}, # mean nir albedo with weak cosz dependency
               'facsf':  {'longName':'facsf','units':'none','ndims':'2','value':0}, # fractional coverage with strong cosz dependency
               'facwf':  {'longName':'facwf','units':'none','ndims':'2','value':0}, # fractional coverage with   weak cosz dependency
               'vfrac':  {'longName':'vfrac','units':'none','ndims':'2','value':tmpvfrac}, # vegetation fraction
               'canopy': {'longName':'canopy','units':'none','ndims':'2','value':0}, # canopy water
               'f10m':   {'longName':'f10m','units':'none','ndims':'2','value':0.997}, # fm at 10 m : Ratio of sigma level 1 wind and 10m wind
               't2m':    {'longName':'t2m','units':'none','ndims':'2','value':tem2m}, # 2 meter temperature
               'q2m':    {'longName':'q2m','units':'none','ndims':'2','value':q2mtmp}, # 2 meter humidity
               'vtype':  {'longName':'vtype','units':'none','ndims':'2','value':tmpvtype}, # vegetation type
               'stype':  {'longName':'stype','units':'none','ndims':'2','value':tmpstype}, # soil type
               'uustar': {'longName':'uustar','units':'none','ndims':'2','value':ust}, # boundary layer parameter (must be nonzero positive)
               'ffmm':   {'longName':'ffmm','units':'none','ndims':'2','value':10}, # fm parameter from PBL scheme
               'ffhh':   {'longName':'ffhh','units':'none','ndims':'2','value':0}, # fh parameter from PBL scheme
               'hice':   {'longName':'hice','units':'none','ndims':'2','value':0}, # sea ice thickness
               'fice':   {'longName':'fice','units':'none','ndims':'2','value':0}, # ice fraction over open water grid
               'tisfc':  {'longName':'tisfc','units':'none','ndims':'2','value':290},# surface temperature over ice fraction
               'tprcp':  {'longName':'tprcp','units':'none','ndims':'2','value':0}, # sfc_fld%tprcp - total precipitation
               'srflag': {'longName':'srflag','units':'none','ndims':'2','value':0},# sfc_fld%srflag - snow/rain flag for precipitation
               'snwdph': {'longName':'snwdph','units':'none','ndims':'2','value':0}, # snow depth water equivalent in mm ; same as snwdph
               'shdmin': {'longName':'shdmin','units':'none','ndims':'2','value':0}, # min fractional coverage of green veg
               'shdmax': {'longName':'shdmax','units':'none','ndims':'2','value':0}, # max fractnl cover of green veg (not used)
               'slope':  {'longName':'slope','units':'none','ndims':'2','value':0}, # sfc slope type for lsm
               'snoalb': {'longName':'snoalb','units':'none','ndims':'2','value':0}, # maximum snow albedo in fraction
               'stc':    {'longName':'stc','units':'none','ndims':'3','value':tmpstc}, # soil temperature (noah)
               'smc':    {'longName':'stc','units':'none','ndims':'3','value':tmpsmc}, # total soil moisture content (noah)
               'slc':    {'longName':'stc','units':'none','ndims':'3','value':tmpslc}, # liquid soil moisture content (noah)
               'tref':   {'longName':'tref','units':'none','ndims':'2','value':290}, # stuff for near surface sea temperature model (NSSTM)
               'z_c':    {'longName':'z_c','units':'none','ndims':'2','value':0.001}, # NSSTM
               'c_0':    {'longName':'c_0','units':'none','ndims':'2','value':0.05}, # NSSTM
               'c_d':    {'longName':'c_d','units':'none','ndims':'2','value':-50}, # NSSTM
               'w_0':    {'longName':'w_0','units':'none','ndims':'2','value':0}, # NSSTM
               'w_d':    {'longName':'w_d','units':'none','ndims':'2','value':0}, # NSSTM
               'xt':     {'longName':'xt','units':'none','ndims':'2','value':0}, # NSSTM
               'xs':     {'longName':'xs','units':'none','ndims':'2','value':0}, # NSSTM
               'xu':     {'longName':'xu','units':'none','ndims':'2','value':0}, # NSSTM
               'xv':     {'longName':'xv','units':'none','ndims':'2','value':0}, # NSSTM
               'xz':     {'longName':'xz','units':'none','ndims':'2','value':30}, # NSSTM
               'zm':     {'longName':'zm','units':'none','ndims':'2','value':0}, # NSSTM
               'xtts':   {'longName':'xtts','units':'none','ndims':'2','value':0}, # NSSTM
               'xzts':   {'longName':'xzts','units':'none','ndims':'2','value':0}, # NSSTM
               'd_conv': {'longName':'d_conv','units':'none','ndims':'2','value':0}, # NSSTM
               'ifd':    {'longName':'ifd','units':'none','ndims':'2','value':0}, # NSSTM
               'dt_cool':{'longName':'dt_cool','units':'none','ndims':'2','value':0.3}, # NSSTM
               'qrain':  {'longName':'qrain','units':'none','ndims':'2','value':0} # NSSTM
               }

oro_vars = {     #  keys     
'slmsk': {'longName':'slmsk','value':tmpslmsk}, # sea/land mask array (sea:0,land:1,sea-ice:2)
'land_frac': {'longName':'land_frac','value':1}, # land  fraction [0:1]
'orog_raw': {'longName':'orog_raw','value':0}, # oro_uf (:)   => null()  !< unfiltered orography
'orog_filt': {'longName':'orog_filt','value':0}, # oro    (:)   => null()  !< orography
'stddev': {'longName':'stddev','value':0}, # hprime (:,1) ; orographic metrics
'convexity': {'longName':'convexity','value':0}, # hprime (:,2) ; orographic metrics
'oa1': {'longName':'oa1','value':0}, # hprime (:,3) ; orographic metrics
'oa2': {'longName':'oa2','value':0}, # hprime (:,4) ; orographic metrics
'oa3': {'longName':'oa3','value':0}, # hprime (:,5) ; orographic metrics
'oa4': {'longName':'oa4','value':0}, # hprime (:,6) ; orographic metrics
'ol1': {'longName':'ol1','value':0}, # hprime (:,7) ; orographic metrics
'ol2': {'longName':'ol2','value':0}, # hprime (:,8) ; orographic metrics
'ol3': {'longName':'ol3','value':0}, # hprime (:,9) ; orographic metrics
'ol4': {'longName':'ol4','value':0}, # hprime (:,10) ; orographic metrics
'theta': {'longName':'theta','value':0}, # hprime (:,11) ; orographic metrics
'gamma': {'longName':'gamma','value':0}, # hprime (:,12) ; orographic metrics
'sigma': {'longName':'sigma','value':0}, # hprime (:,13) ; orographic metrics
'elvmax': {'longName':'elvmax','value':0} # hprime (:,14) ; orographic metrics
}

ncfile = Dataset(sfcfilename,mode='w',clobber=False,format='NETCDF4_CLASSIC') 


x_dim = ncfile.createDimension('xaxis_1', nx)     # latitude axis
y_dim = ncfile.createDimension('yaxis_1', ny)    # longitude axis
z_dim = ncfile.createDimension('zaxis_1', nz)    # longitude axis
time_dim = ncfile.createDimension('Time', 1) # unlimited axis (can be appended to).

xdim = ncfile.createVariable('xaxis_1', np.float32, ('xaxis_1',))
xdim.units = 'none'
xdim.longName = 'xaxis_1'
xdim.cartesian_axis = 'X'
xdim[:] = 1 + np.arange(nx)

ydim = ncfile.createVariable('yaxis_1', np.float32, ('yaxis_1',))
ydim.units = 'none'
ydim.longName = 'yaxis_1'
ydim.cartesian_axis = 'Y'
ydim[:] = 1 + np.arange(ny)

zdim = ncfile.createVariable('zaxis_1', np.float32, ('zaxis_1',))
zdim.units = 'none'
zdim.longName = 'zaxis_1'
zdim.cartesian_axis = 'Z'
zdim[:] = 1 + np.arange(nz)

timedim = ncfile.createVariable('Time', np.float32, ('Time',))
timedim.units = 'Time level'
timedim.longName = 'Time'
timedim.cartesian_axis = 'T'
timedim[:] = 1

geolat = ncfile.createVariable('geolat', np.float64, ('yaxis_1','xaxis_1'),zlib=True)
geolat.units = 'degrees_north'
geolat.longName = 'Latitude'
geolat[:] = lat1

geolon = ncfile.createVariable('geolon', np.float64, ('yaxis_1','xaxis_1'),zlib=True)
geolon.units = 'degrees_east'
geolon.longName = 'Longitude'
geolon[:] = lon1

for k, key in enumerate(sfc_vars):
    if sfc_vars[key]["ndims"] == '2':
       var = ncfile.createVariable(key, np.float64, ('yaxis_1','xaxis_1'),zlib=True)
    else:
       var = ncfile.createVariable(key, np.float64, ('zaxis_1','yaxis_1','xaxis_1'),zlib=True)
    #print('key = ',key)
    print('sfcvar = ',sfc_vars[key])
    #print('sfcvar = ',sfc_vars[key]["longName"])
    var.longName = key
    var.units = sfc_vars[key]["units"]
    var.coordinates = "geolon geolat" 
    var[:] = sfc_vars[key].get("value")


# first print the Dataset object to see what we've got
print(ncfile)
# close the Dataset.
ncfile.close(); print('Sfc Dataset is closed!')


# Ideal oro_data.nc
ncfile = Dataset(orofilename,mode='w',clobber=False,format='NETCDF4_CLASSIC') 

lat = ncfile.createDimension('lat', nx)     # latitude axis
lon = ncfile.createDimension('lon', ny)    # longitude axis

geolat = ncfile.createVariable('geolat', np.float32, ('lat','lon'),zlib=True)
geolat.units = 'degrees_north'
geolat.long_name = 'Latitude'
geolat[:] = lat1

geolon = ncfile.createVariable('geolon', np.float32, ('lat','lon'),zlib=True)
geolon.units = 'degrees_east'
geolon.long_name = 'Longitude'
geolon[:] = lon1


for k, key in enumerate(oro_vars):
    var = ncfile.createVariable(key, np.float32, ('lat','lon'),zlib=True)
    #print('key = ',key)
    print('orovar = ',oro_vars[key])
    #print('sfcvar = ',sfc_vars[key]["longName"])
    var.coordinates = "geolon geolat"
    var[:] = oro_vars[key].get("value")


# first print the Dataset object to see what we've got
print(ncfile)
# close the Dataset.
ncfile.close(); print('Orog Dataset is closed!')

