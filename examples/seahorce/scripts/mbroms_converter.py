import os
import sys
import netCDF4
import numpy as np
from optparse import OptionParser

# sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append("/opt/moab/seahorce/lib/python3.9/site-packages")

from pymoab import core
from pymoab import types
from pymoab.rng import Range

# print("== Gathering information.  (Invoke with --help for more details. All arguments are optional)")
# parser = OptionParser()
# parser.description = "This script takes a scrip file and generates a MOAB h5m file."

# parser.add_option("-m", "--moab", dest="moabFile", help="MOAB ROMS grid file name input with mesh and field data.", default="roms_3d_projected.h5m", metavar="FILENAME")
# parser.add_option("-r", "--roms", dest="romsFile", help="ROMS grid file output used for restart.", default="roms_his.nc", metavar="FILENAME")

# for option in parser.option_list:
#   if option.default != ("NO", "DEFAULT"):
#       option.help += (" " if option.help else "") + "[default: %default]"
                
# options, args = parser.parse_args()

# if not options.moabFile:
#   sys.exit('Error: MOAB output grid filename is required.  Specify with -m command line argument.')
# if not options.romsFile:
#   sys.exit('Error: ROMS output grid filename is required.  Specify with -r command line argument.')

# moabInputFile = options.moabFile
# romsOutputFile = options.romsFile

moabInputFile = "roms_2d_projected.h5m"
romsOutputFile = "roms_his.nc"

# start a MOAB instance
mb1 = core.Core()
mb2 = core.Core()

# load the file
twod_set = mb1.create_meshset()
mb1.load_file("roms_2d_projected.h5m", twod_set)
threed_set = mb2.create_meshset()
mb2.load_file("roms_3d_projected.h5m", threed_set)

# query the root set for all vertices
verts2d = mb1.get_entities_by_type(twod_set, types.MBVERTEX, recur = True)
print("Found " + str(verts2d.size()) + " vertices in this model.")

ents2d = mb1.get_entities_by_dimension(twod_set, 2)
print("Found " + str(ents2d.size()) + " 2D faces in this model.")

ents3d = mb2.get_entities_by_dimension(threed_set, 3)
print("Found " + str(ents3d.size()) + " 3D elements in this model.")

try: ncfile.close()  # just to be safe, make sure dataset is not already open.
except: pass

ncfile = netCDF4.Dataset(romsOutputFile,mode='w',format='NETCDF4_CLASSIC')
print(ncfile)

# 2D elements: 154 * 142 rho points
xi_rho_init = 154
eta_rho_init = 142
s_rho_init = int(ents3d.size()/ents2d.size())

assert(xi_rho_init*eta_rho_init == ents2d.size())
# assert(s_rho_init == ents3d.size()/ents2d.size())

salinity_tag = mb2.tag_get_handle("Salinity",1,types.MB_TYPE_DOUBLE,types.MB_TAG_DENSE)
# print("Salinity: ", salinity_tag)
# make sure we can still get data for the write tag
#salinity = mb2.tag_get_data(salinity_tag, ents3d).reshape(xi_rho_init, eta_rho_init, s_rho_init)
salinity = mb2.tag_get_data(salinity_tag, ents3d).reshape(s_rho_init, xi_rho_init, eta_rho_init)
# salinity = mb2.tag_get_data(salinity_tag, ents3d)
# print("Salinity: ", salinity_tag, salinity.shape)
# sal3d = salinity[:,0].reshape(xi_rho_init, eta_rho_init, s_rho_init)

temperature_tag = mb2.tag_get_handle("Temperature",1,types.MB_TYPE_DOUBLE,types.MB_TAG_DENSE)
#temperature = mb2.tag_get_data(temperature_tag, ents3d).reshape(xi_rho_init, eta_rho_init, s_rho_init)
temperature = mb2.tag_get_data(temperature_tag, ents3d).reshape(s_rho_init, xi_rho_init, eta_rho_init)

xi_rho = ncfile.createDimension('xi_rho', xi_rho_init) # longitude axis
eta_rho = ncfile.createDimension('eta_rho', eta_rho_init) # longitude axis
s_rho = ncfile.createDimension('s_rho', s_rho_init) # latitude axis
time_dim = ncfile.createDimension('ocean_time', None) # unlimited axis (can be appended to).
for dim in ncfile.dimensions.items():
  print(dim)

# attributes
ncfile.title='Wind-Driven Upwelling/Downwelling for the Niskine test case'
ncfile.type="ROMS/TOMS restart file"

print(ncfile.title)
print(ncfile.type)
print(ncfile)

## variables
# free surface height: double zeta(ocean_time, eta_rho, xi_rho) ;
# zeta = ncfile.createVariable('zeta', np.float64, ('ocean_time','eta_rho', 'xi_rho'))
# zeta.long_name = "free-surface" ;
# zeta.units = "meter" ;
# zeta.time = "ocean_time" ;
# zeta.grid = "grid" ;
# zeta.location = "face" ;
# zeta.coordinates = "x_rho y_rho ocean_time" ;
# zeta.field = "free-surface, scalar, series" ;
# print(zeta)

temp = ncfile.createVariable('temp',np.float64,('ocean_time','s_rho','eta_rho','xi_rho')) # note: unlimited dimension is leftmost
temp.units = 'Celsius' # degrees Celsius
temp.long_name = 'potential temperature' # this is a CF standard name
temp.time = "ocean_time" ;
temp.grid = "grid" ;
temp.location = "face" ;
temp.coordinates = "x_rho y_rho s_rho ocean_time" ;
temp.field = "temperature, scalar, series" ;
print(temp)
#print("Found shape: ", temperature.shape, np.einsum('ijk->kij', temperature).shape)
temp[0, :, :, :] = temperature
#temp[0, :, :, :] = np.einsum('ijk->kij', temperature)
print(temp)

salt = ncfile.createVariable('salt',np.float64,('ocean_time','s_rho','eta_rho','xi_rho')) # note: unlimited dimension is leftmost
salt.long_name = 'salinity' # this is a CF standard name
salt.time = "ocean_time" ;
salt.grid = "grid" ;
salt.location = "face" ;
salt.coordinates = "x_rho y_rho s_rho ocean_time" ;
salt.field = "salinity, scalar, series" ;
salt[0, :, :, :] = salinity
#salt[0, :, :, :] = np.einsum('ijk->kij', salinity)
print(salt)

# custom attributes
print("-- Some pre-defined attributes for variable temp:")
print("temp.dimensions:", temp.dimensions)
print("temp.shape:", temp.shape)
print("temp.dtype:", temp.dtype)
print("temp.ndim:", temp.ndim)

## write
# nlats = len(lat_dim); nlons = len(lon_dim); ntimes = 3
# lat[:] = -90. + (180./nlats)*np.arange(nlats) # south pole to north pole
# lon[:] = (180./nlats)*np.arange(nlons) # Greenwich meridian eastward
# data_arr = np.random.uniform(low=280,high=330,size=(ntimes,nlats,nlons))
# temp[:,:,:] = data_arr # Appends data along unlimited dimension
# print("-- Wrote data, temp.shape is now ", temp.shape)
# print("-- Min/Max values:", temp[:,:,:].min(), temp[:,:,:].max())

# data_slice = np.random.uniform(low=280,high=330,size=(nlats,nlons))
# temp[3,:,:] = data_slice
# print("-- Wrote more data, temp.shape is now ", temp.shape)

# print(time)
# times_arr = time[:]
# print(type(times_arr),times_arr)


# first print the Dataset object to see what we've got
print(ncfile)
# close the Dataset.
ncfile.close(); print('Dataset is closed!')

