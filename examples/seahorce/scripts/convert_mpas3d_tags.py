import os
import sys

# sys.path.append("/opt/moab/seahorce/lib/python3.9/site-packages")
sys.path.append("/Users/mahadevan/Library/Python/3.9/lib/python/site-packages")

import pymoab
import numpy as np
from pymoab import core,types
import netCDF4

## Prepricess h5m file:
# ../../build/tools/mbconvert -O "VARIABLE=" -O "NO_EDGES" -O "NO_MIXED_ELEMENTS"  mpas_grid.nc mpas_grid_raw.h5m

## List of tags that we care about
## Only z: refBottomDepth, refZMid, vertCoordMovementWeights
## Only edges: dvEdge, edgeMask(z)
## Cells/z: restingThickness
taglist_clean = ["nCells", "nEdges", "nVertices", "nVertLevels", "maxEdges", "maxEdges2", "vertexDegree", "Time", "TWO", "bedElevation", "dvEdge",  "bed_elevation",  "edgeMask", "bottomDepthObserved", "boundaryLayerDepth0", "fCell", "refBottomDepth",  "refZMid", "restingThickness", "salinity0", "temperature0", "vertCoordMovementWeights", "vertexDegree", "layerThickness0"]
taglist_d = ["bottomDepth"] # ["bed_elevation", "bottomDepth", "bottomDepthObserved", "fCell"]
taglist_i = ["maxLevelCell", "minLevelCell"]
taglist_z = ["refBottomDepth"]
taglist_4d = ["temperature", "salinity", "layerThickness"]
taglist_42d = ["temperature", "salinity"]
forcing_data = ["seaIcePressure", "atmosphericPressure", "windStressZonal", "windStressMeridional", "evaporationFlux", "rainFlux", "temperaturePistonVelocity", "salinityPistonVelocity", "temperatureSurfaceRestoringValue", "salinitySurfaceRestoringValue"]
# hdf5_filename = "mpas_grid_raw.h5m"
hdf5_filename = "mpaso_30to10km_standalone_rerun/mpas_mesh_test.h5m"
ncforce_filename = "mpaso_30to10km_standalone_rerun/forcing_data.nc"
nc_filename = "mpaso_30to10km_standalone_rerun/mpas_init_data.nc"

# start a MOAB instance
mb = core.Core()

# load the file
mb.load_file(hdf5_filename)

# get the root set of the MOAB instance
root_set = mb.get_root_set()

gidTag = mb.tag_get_handle("GLOBAL_ID")

# query the root set for all polygons
polys = mb.get_entities_by_type(root_set, types.MBPOLYGON, recur = True)
print("Found " + str(polys.size()) + " polygons in this model.")

# similar query for vertices
verts = mb.get_entities_by_type(root_set, types.MBVERTEX, recur = True)
print("Found " + str(verts.size()) + " vertices in this model.\n")

gids = mb.tag_get_data(gidTag,polys).reshape((polys.size())) - 1

taglist_clean.extend(taglist_d)
taglist_clean.extend(taglist_i)
taglist_clean.extend(taglist_z)
taglist_clean.extend(taglist_4d)
taglist_clean.extend(taglist_42d)
print("Taglist: ", taglist_clean)
for etag in taglist_clean:
    print("Deleting", etag)
    try:
        thandle = mb.tag_get_handle(etag)
        mb.tag_delete(thandle)
    except Exception as e:
        pass

tdata = np.zeros((polys.size()))

# get the tag data
ncf = netCDF4.Dataset(nc_filename, 'r')
ncforce = netCDF4.Dataset(ncforce_filename, 'r')
nVertLevels = ncf.dimensions["nVertLevels"].size

print("Global IDs: ", gids[:10])

for dtag in taglist_d:
    print("\nAnalyzing", dtag)

    ncvar = ncf.variables[dtag]
    print(ncvar)

    # get the actual data out of the variable
    tdata[gids[:]] = ncvar[:]

    # print(tdata[:10])

    print("Setting", dtag, "tag data")
    thandle = mb.tag_get_handle(dtag,1,types.MB_TYPE_DOUBLE,types.MB_TAG_DENSE,True)
    mb.tag_set_data(thandle,polys,tdata)

idata = np.zeros((polys.size()), dtype=np.int32)

for itag in taglist_i:
    print("\nAnalyzing", itag)

    ncvar = ncf.variables[itag]
    print(ncvar)

    # get the actual data out of the variable
    idata[gids[:]] = ncvar[:]

    print("Setting", itag, "tag data")
    thandle = mb.tag_get_handle(itag,1,types.MB_TYPE_INTEGER,types.MB_TAG_DENSE,True)
    mb.tag_set_data(thandle,polys,idata)

zdata = np.zeros((nVertLevels))
for ztag in taglist_z:
    print("\nAnalyzing", ztag)

    ncvar = ncf.variables[ztag]
    print(ncvar)

    # get the actual data out of the variable
    zdata[:] = ncvar[:]

    print("Setting", ztag, "tag data")
    thandle = mb.tag_get_handle(ztag,nVertLevels,types.MB_TYPE_DOUBLE,types.MB_TAG_SPARSE,True)
    mb.tag_set_data(thandle,root_set,zdata)

# tdata3d = np.zeros((polys.size(), nVertLevels))
tdata3d = np.zeros((polys.size() * nVertLevels))

for fdtag in taglist_4d:
    print("\nAnalyzing", fdtag)

    ncvar = ncf.variables[fdtag]
    print(ncvar)

    # get the actual data out of the variable
    tdata3d[:] = ncvar[0, :, :].flatten()

    # print(tdata3d[300*60:301*60])
    # print(tdata3d[301*60:302*60])

    print("Setting", fdtag, "tag data")
    thandle = mb.tag_get_handle(fdtag+"_3d",nVertLevels,types.MB_TYPE_DOUBLE,types.MB_TAG_DENSE,True)
    mb.tag_set_data(thandle,polys,tdata3d)

for fdtag in taglist_42d:
    print("\nAnalyzing", fdtag)

    ncvar = ncf.variables[fdtag]
    print(ncvar)

    # get the actual data out of the variable
    tdata[:] = ncvar[0, :, 0].flatten()

    print("Setting", fdtag, "tag data")
    thandle = mb.tag_get_handle(fdtag+"_2d",1,types.MB_TYPE_DOUBLE,types.MB_TAG_DENSE,True)
    mb.tag_set_data(thandle,polys,tdata)

for dtag in forcing_data:
    print("\nAnalyzing", dtag)

    ncvar = ncforce.variables[dtag]
    print(ncvar, tdata.shape)

    # get the actual data out of the variable
    tdata[gids[:]] = ncvar[:]

    # print(tdata[:10])
    print("Setting", dtag, "tag data")
    thandle = mb.tag_get_handle(dtag,1,types.MB_TYPE_DOUBLE,types.MB_TAG_DENSE,True)
    mb.tag_set_data(thandle,polys,tdata)

ncf.close()
ncforce.close()
mb.write_file("mpas_grid_with_fields.h5m")
