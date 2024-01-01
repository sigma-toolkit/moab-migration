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

niskine = False
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

if niskine:
    titlestr='Wind-Driven Upwelling/Downwelling for the Niskine test case'
    # 2D elements: 154 * 142 rho points
    xi_rho_init = 154
    eta_rho_init = 142
    xi_u_init = 153
    eta_u_init = 142
    xi_v_init = 154
    eta_v_init = 141
    s_rho_init = int(ents3d.size()/ents2d.size())
else:
    titlestr='Wind-Driven Upwelling/Downwelling for the North-Atlantic test case'
    # 2D elements: 154 * 142 rho points
    xi_rho_init = 413
    eta_rho_init = 147
    xi_u_init = 412
    eta_u_init = 147
    xi_v_init = 413
    eta_v_init = 146
    #s_rho_init = 100
    s_rho_init = int(ents3d.size()/ents2d.size())

assert(xi_rho_init*eta_rho_init == ents2d.size())
# assert(s_rho_init == ents3d.size()/ents2d.size())

bathmetry_tag = mb1.tag_get_handle("Bathymetry",1,types.MB_TYPE_DOUBLE,types.MB_TAG_DENSE)
bathymetry = mb1.tag_get_data(bathmetry_tag, ents2d).reshape(eta_rho_init, xi_rho_init)

ssh_tag = mb1.tag_get_handle("SSH",1,types.MB_TYPE_DOUBLE,types.MB_TAG_DENSE)
ssh = mb1.tag_get_data(ssh_tag, ents2d).reshape(eta_rho_init, xi_rho_init)

salinity_tag = mb2.tag_get_handle("ROMS_Salinity",1,types.MB_TYPE_DOUBLE,types.MB_TAG_DENSE)
salinity = mb2.tag_get_data(salinity_tag, ents3d).reshape(s_rho_init, xi_rho_init, eta_rho_init)[::-1,:,:]
# salinity = mb2.tag_get_data(salinity_tag, ents3d).reshape(s_rho_init, eta_rho_init, xi_rho_init)[::-1,:,:]
print("Salinity: ", salinity.shape)

temperature_tag = mb2.tag_get_handle("ROMS_Temperature",1,types.MB_TYPE_DOUBLE,types.MB_TAG_DENSE)
temperature = mb2.tag_get_data(temperature_tag, ents3d).reshape(s_rho_init, xi_rho_init, eta_rho_init)[::-1,:,:]

velocityx_tag = mb2.tag_get_handle("ROMS_VelZonal",1,types.MB_TYPE_DOUBLE,types.MB_TAG_DENSE)
velocityx = mb2.tag_get_data(velocityx_tag, ents3d).reshape(s_rho_init, eta_rho_init, xi_rho_init)[::-1,:,:]
velocityxavg = 0.5 * ( velocityx[:,:,:-1] + velocityx[:,:,1:] )
# velocityxavg = ( velocityx[:,:,:-1] )
print("X velocity shape: ", velocityx.shape, velocityxavg.shape)

velocityy_tag = mb2.tag_get_handle("ROMS_VelMeridional",1,types.MB_TYPE_DOUBLE,types.MB_TAG_DENSE)
velocityy = mb2.tag_get_data(velocityy_tag, ents3d).reshape(s_rho_init, eta_rho_init, xi_rho_init)[::-1,:,:]
velocityyavg = 0.5 * ( velocityy[:,:-1,:] + velocityy[:,1:,:] )
# velocityyavg = ( velocityy[:,:-1,:] )
print("Y velocity shape: ", velocityy.shape, velocityyavg.shape)

xi_rho = ncfile.createDimension('xi_rho', xi_rho_init) # longitude axis
eta_rho = ncfile.createDimension('eta_rho', eta_rho_init) # longitude axis
s_rho = ncfile.createDimension('s_rho', s_rho_init) # vertical axis
xi_u = ncfile.createDimension('xi_u', xi_u_init) # longitude axis
eta_u = ncfile.createDimension('eta_u', eta_u_init) # longitude axis
xi_v = ncfile.createDimension('xi_v', xi_v_init) # longitude axis
eta_v = ncfile.createDimension('eta_v', eta_v_init) # longitude axis
time_dim = ncfile.createDimension('ocean_time', None) # unlimited axis (can be appended to).
for dim in ncfile.dimensions.items():
  print(dim)

# attributes
ncfile.title=titlestr
ncfile.type="ROMS/TOMS restart file"

print(ncfile.title)
print(ncfile.type)
# print(ncfile)

## variables
theta_s = ncfile.createVariable('theta_s',np.float64,())
theta_s.long_name = "S-coordinate surface control parameter"
theta_s = 5.0

theta_b = ncfile.createVariable('theta_b',np.float64,())
theta_b.long_name = "S-coordinate bottom control parameter"
theta_b = 0.5

# Tcline = ncfile.createVariable('Tcline',np.float64,())
# Tcline.long_name = "S-coordinate surface/bottom layer width"
# Tcline.units = "meter"
# Tcline = 0.5

hc = ncfile.createVariable('hc',np.float64,())
hc.long_name = "S-coordinate parameter, critical depth"
hc.units = "meter"
hc = 100.0

        # double s_rho(s_rho) ;
        #         s_rho:long_name = "S-coordinate at RHO-points" ;
        #         s_rho:valid_min = -1. ;
        #         s_rho:valid_max = 0. ;
        #         s_rho:positive = "up" ;
        #         s_rho:standard_name = "ocean_s_coordinate_g2" ;
        #         s_rho:formula_terms = "s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc" ;
        #         s_rho:field = "s_rho, scalar" ;
        # double s_w(s_w) ;
        #         s_w:long_name = "S-coordinate at W-points" ;
        #         s_w:valid_min = -1. ;
        #         s_w:valid_max = 0. ;
        #         s_w:positive = "up" ;
        #         s_w:standard_name = "ocean_s_coordinate_g2" ;
        #         s_w:formula_terms = "s: s_w C: Cs_w eta: zeta depth: h depth_c: hc" ;
        #         s_w:field = "s_w, scalar" ;
        # double Cs_r(s_rho) ;
        #         Cs_r:long_name = "S-coordinate stretching curves at RHO-points" ;
        #         Cs_r:valid_min = -1. ;
        #         Cs_r:valid_max = 0. ;
        #         Cs_r:field = "Cs_r, scalar" ;
        # double Cs_w(s_w) ;
        #         Cs_w:long_name = "S-coordinate stretching curves at W-points" ;
        #         Cs_w:valid_min = -1. ;
        #         Cs_w:valid_max = 0. ;
        #         Cs_w:field = "Cs_w, scalar" ;


# Bathymetry field: double h(eta_rho, xi_rho) ;
bath = ncfile.createVariable('h',np.float64,('eta_rho','xi_rho'))
bath.long_name = "bathymetry at RHO-points" ;
bath.units = "meter" ;
bath.grid = "grid" ;
bath.location = "face" ;
bath.coordinates = "lon_rho lat_rho" ;
bath.field = "bath, scalar" ;
bath[:,:] = bathymetry
print(bath)

# free surface height: double zeta(ocean_time, eta_rho, xi_rho) ;
zeta = ncfile.createVariable('zeta', np.float64, ('ocean_time','eta_rho', 'xi_rho'))
zeta.long_name = "free-surface" ;
zeta.units = "meter" ;
zeta.time = "ocean_time" ;
zeta.grid = "grid" ;
zeta.location = "face" ;
zeta.coordinates = "x_rho y_rho ocean_time" ;
zeta.field = "free-surface, scalar, series" ;
zeta[0,:,:] = ssh
print(zeta)

fieldvars = {
             'salinityPistonVelocity': "m s^-1",
             'temperaturePistonVelocity': "m s^-1",
             'salinitySurfaceRestoringValue': "C",
             'temperatureSurfaceRestoringValue': "C",
             'windStressMeridional': "N m^-2",
             'windStressZonal': "N m^-2"
             }

fielddesc = {
             'salinityPistonVelocity': "A non-negative field controlling the rate at which salinity is restored to salinitySurfaceRestoringValue",
             'temperaturePistonVelocity': "A non-negative field controlling the rate at which temperature is restored to temperatureSurfaceRestoringValue",
             'salinitySurfaceRestoringValue': "Salinity is restored toward this field at a rate controlled by salinityPistonVelocity.",
             'temperatureSurfaceRestoringValue': "Temperature is restored toward this field at a rate controlled by temperaturePistonVelocity.",
             'windStressMeridional': "Meridional (northward) component of wind stress at cell centers from coupler. Positive northward.",
             'windStressZonal': "Zonal (eastward) component of wind stress at cell centers from coupler. Positive eastward."
             }

# Bathymetry field: double h(eta_rho, xi_rho) ;
for fieldname in fieldvars:

    field_tag = mb1.tag_get_handle(fieldname,1,types.MB_TYPE_DOUBLE,types.MB_TAG_DENSE)
    if field_tag:
      fieldata = mb1.tag_get_data(field_tag, ents2d).reshape(eta_rho_init, xi_rho_init)

      field = ncfile.createVariable(fieldname,np.float64,('eta_rho','xi_rho'))
      field.long_name = fielddesc[fieldname] ;
      field.units = fieldvars[fieldname] ;
      field.grid = "grid" ;
      field.location = "face" ;
      field.coordinates = "lon_rho lat_rho" ;
      field.field = fieldname + ", scalar" ;
      field[:,:] = fieldata
      print(field)

# Temperature scalar field
temp = ncfile.createVariable('temp',np.float64,('ocean_time','s_rho','eta_rho','xi_rho')) # note: unlimited dimension is leftmost
temp.units = 'Celsius' # degrees Celsius
temp.long_name = 'potential temperature' # this is a CF standard name
temp.time = "ocean_time" ;
temp.grid = "grid" ;
temp.location = "face" ;
temp.coordinates = "x_rho y_rho s_rho ocean_time" ;
temp.field = "temperature, scalar, series" ;
print("Temperature shape: ", temperature.shape)
#print("Found shape: ", temperature.shape, np.einsum('ijk->kij', temperature).shape)
temp[0, :, :, :] = temperature
#temp[0, :, :, :] = np.einsum('ijk->kij', temperature)
print(temp)

# Salinity scalar field
salt = ncfile.createVariable('salt',np.float64,('ocean_time','s_rho','eta_rho','xi_rho')) # note: unlimited dimension is leftmost
salt.long_name = 'salinity' # this is a CF standard name
salt.time = "ocean_time" ;
salt.grid = "grid" ;
salt.location = "face" ;
salt.coordinates = "x_rho y_rho s_rho ocean_time" ;
salt.field = "salinity, scalar, series" ;
print("Salinity shape: ", salinity.shape)
salt[0, :, :, :] = salinity
#salt[0, :, :, :] = np.einsum('ijk->kij', salinity)
print(salt)

# Velocity fields
# float u(ocean_time, s_rho, eta_u, xi_u) ;
velu = ncfile.createVariable('u',np.float64,('ocean_time','s_rho','eta_u','xi_u')) # note: unlimited dimension is leftmost
# velu = ncfile.createVariable('u',np.float64,('ocean_time','s_rho','eta_rho','xi_rho')) # note: unlimited dimension is leftmost
velu.long_name = "u-momentum component" ;
velu.units = "meter second-1" ;
velu.time = "ocean_time" ;
velu.grid = "grid" ;
velu.location = "edge1" ;
velu.coordinates = "lon_u lat_u s_rho ocean_time" ;
velu.field = "u-velocity, scalar, series" ;
print("X velocity shape: ", velocityxavg.shape, np.delete(velocityx, 0, 1).shape)
# velu[0, :, :, :] = velocityx
# velu[0, :, :, :] = np.delete(velocityx, 0, 1)
velu[0, :, :, :] = velocityxavg
#velu[0, :, :, :] = np.einsum('ijk->ikj', velocityxavg)
print(velu)

# float v(ocean_time, s_rho, eta_v, xi_v) ;
velv = ncfile.createVariable('v',np.float64,('ocean_time','s_rho','eta_v','xi_v')) # note: unlimited dimension is leftmost
# velv = ncfile.createVariable('v',np.float64,('ocean_time','s_rho','eta_rho','xi_rho')) # note: unlimited dimension is leftmost
velv.long_name = "v-momentum component" ;
velv.units = "meter second-1" ;
velv.time = "ocean_time" ;
velv.grid = "grid" ;
velv.location = "edge2" ;
velv.coordinates = "lon_v lat_v s_rho ocean_time" ;
velv.field = "v-velocity, scalar, series" ;
print("Y velocity shape: ", velocityyavg.shape, np.delete(velocityy, 0, 2).shape)
print(velv)
# velv[0, :, :, :] = velocityy
# velv[0, :, :, :] = np.delete(velocityy, 0, 2)
velv[0, :, :, :] = velocityyavg
# velv[0, :, :, :] = np.einsum('ijk->ikj', velocityyavg)
print(velv)

# custom attributes
# print("-- Some pre-defined attributes for variable temp:")
# print("temp.dimensions:", temp.dimensions)
# print("temp.shape:", temp.shape)
# print("temp.dtype:", temp.dtype)
# print("temp.ndim:", temp.ndim)

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
print("\n\n", ncfile)
# close the Dataset.
ncfile.close(); print('Dataset is closed!')

