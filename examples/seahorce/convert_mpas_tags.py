import numpy as np
from pymoab import core,types
import netCDF4

## Prepricess h5m file:
# ../../build/tools/mbconvert -O "VARIABLE=nCells" -O "NO_EDGES" -O "NO_MIXED_ELEMENTS"  mpas_grid.nc mpas_grid_raw.h5m

## List of tags that we care about
## Only z: refBottomDepth, refZMid, vertCoordMovementWeights
## Only edges: dvEdge, edgeMask(z)
## Cells/z: restingThickness
taglist_clean = ["nCells", "nEdges", "nVertices", "nVertLevels", "maxEdges", "maxEdges2", "vertexDegree", "Time", "TWO"]
taglist_d = ["bottomDepth"] # ["bed_elevation", "bottomDepth", "bottomDepthObserved", "fCell"]
taglist_i = [] # ["maxLevelCell", "minLevelCell"]
taglist_z = ["refBottomDepth"]
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

# get the tag data
ncf = netCDF4.Dataset(nc_filename, 'r')
nVertLevels = ncf.dimensions["nVertLevels"].size
for dtag in taglist_d:
    print("\nAnalyzing", dtag)

    ncvar = ncf.variables[dtag]
    print(ncvar)

    # get the actual data out of the variable
    tdata[gids[:]] = ncvar[:]

    print("Setting", dtag, "tag data")
    thandle = mb.tag_get_handle(dtag,1,types.MB_TYPE_DOUBLE,types.MB_TAG_DENSE,True)
    mb.tag_set_data(thandle,polys,tdata)

for itag in taglist_i:
    print("\nAnalyzing", itag)

    ncvar = ncf.variables[itag]
    print(ncvar)

    # get the actual data out of the variable
    tdata[gids[:]] = ncvar[:]

    print("Setting", itag, "tag data")
    thandle = mb.tag_get_handle(itag,1,types.MB_TYPE_INTEGER,types.MB_TAG_DENSE,True)
    mb.tag_set_data(thandle,polys,tdata)

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


ncf.close()
mb.write_file("mpas_grid.h5m")

