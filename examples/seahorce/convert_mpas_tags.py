## Workflow:
# Convert: mbconvert mpasgrid.nc mpas_grid_raw.h5m
# Preprocess: python convert_mpas_tags.py

import numpy as np
from pymoab import core,types
import netCDF4

## Prepricess h5m file:
# ../../build/tools/mbconvert -O "VARIABLE=nCells" -O "NO_EDGES" -O "NO_MIXED_ELEMENTS"  mpas_grid.nc mpas_grid_raw.h5m

## List of tags that we care about
## Only z: refBottomDepth, refZMid, vertCoordMovementWeights
## Only edges: dvEdge, edgeMask(z)
## Cells/z: restingThickness
taglist_clean = ["nCells", "nEdges", "nVertices", "nVertLevels", "maxEdges", "maxEdges2", "vertexDegree", "Time", "TWO", "bed_elevation", "bottomDepth", "bottomDepthObserved", "boundaryLayerDepth0", "dvEdge", "edgeMask", "fCell", "layerThickness0", "maxLevelCell", "minLevelCell", "normalVelocity0", "refBottomDepth", "refZMid", "restingThickness", "salinity0", "temperature0", "vertCoordMovementWeights"]
taglist_d = ["bottomDepth"] # ["bed_elevation", "bottomDepth", "bottomDepthObserved", "fCell"]
taglist_i = ["maxLevelCell", "minLevelCell"]
taglist_z = ["refBottomDepth"]
taglist_4d = ["temperature", "salinity", "layerThickness"]
taglist_42d = ["temperature", "salinity"]
hdf5_filename = "mpas_grid_raw.h5m"
nc_filename = "mpas_grid.nc"

# start a MOAB instance
mb = core.Core()

# load the file
mb.load_file(hdf5_filename)

# get the root set of the MOAB instance
root_set = mb.get_root_set()

gidTag = mb.tag_get_handle("GLOBAL_ID")

# query the root set for all triangles
polys = mb.get_entities_by_type(root_set, types.MBPOLYGON, recur = True)
print("Found " + str(polys.size()) + " polygons in this model.")

# similar query for vertices
verts = mb.get_entities_by_type(root_set, types.MBVERTEX, recur = True)
print("Found " + str(verts.size()) + " vertices in this model.\n")

gids = mb.tag_get_data(gidTag,polys)

gids = gids.reshape((polys.size())) - 1

for etag in taglist_clean:
    print("Deleting", etag)
    thandle = mb.tag_get_handle(etag)
    mb.tag_delete(thandle)

tdata = np.zeros((polys.size()))
idata = np.zeros((polys.size()), dtype=np.int32)

# get the tag data
ncf = netCDF4.Dataset(nc_filename, 'r')
nVertLevels = ncf.dimensions["nVertLevels"].size

if True:
    for dtag in taglist_d:
        print("\nAnalyzing", dtag)

        ncvar = ncf.variables[dtag]
        print(ncvar)

        # get the actual data out of the variable
        # tdata[gids[:]] = ncvar[:]
        tdata[:] = ncvar[gids[:]]

        # print(tdata[:50])

        print("Setting", dtag, "tag data")
        thandle = mb.tag_get_handle(dtag,1,types.MB_TYPE_DOUBLE,types.MB_TAG_DENSE,True)
        mb.tag_set_data(thandle,polys,tdata)

    for itag in taglist_i:
        print("\nAnalyzing", itag)

        ncvar = ncf.variables[itag]
        print(ncvar)

        # get the actual data out of the variable
        idata[:] = ncvar[gids[:]]

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
        tdata3d[:] = ncvar[0, gids[:], :].flatten()

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
        tdata[:] = ncvar[0, gids[:], 0].flatten()

        print("Setting", fdtag, "tag data")
        thandle = mb.tag_get_handle(fdtag,1,types.MB_TYPE_DOUBLE,types.MB_TAG_DENSE,True)
        mb.tag_set_data(thandle,polys,tdata)


## verify
if False:
    ncvar = ncf.variables["layerThickness"]
    print(ncvar)

    # get the actual data out of the variable
    tdata3d_a = ncvar[0, gids[1], :]
    tdata3d_b = ncvar[0, gids[2], :]
    tdata3d_c = ncvar[0, gids[3], :]

    print(tdata3d_a[:]-tdata3d_b[:])
    print(tdata3d_b[:]-tdata3d_c[:])

ncf.close()
mb.write_file("mpas_grid.h5m")



