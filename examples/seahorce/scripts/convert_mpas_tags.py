## Workflow:
# Convert: mbconvert mpasgrid.nc mpas_grid_raw.h5m
# Preprocess: python convert_mpas_tags.py
import sys, getopt
import numpy as np
from pymoab import core, types
import netCDF4

## Prepricess h5m file:
# ../../build/tools/mbconvert -O "VARIABLE=nCells" -O "NO_EDGES" -O "NO_MIXED_ELEMENTS"  mpas_grid.nc mpas_grid_raw.h5m

## List of tags that we care about
## Only z: refBottomDepth, refZMid, vertCoordMovementWeights
## Only edges: dvEdge, edgeMask(z)
## Cells/z: restingThickness
taglist_clean = [
    "nCells",
    "nEdges",
    "nVertices",
    "nVertLevels",
    "maxEdges",
    "maxEdges2",
    "vertexDegree",
    "Time",
    "TWO",
    #"bed_elevation",
    #"bottomDepth",
    #"bottomDepthObserved",
    #"boundaryLayerDepth0",
    #"dvEdge",
    #"edgeMask",
    #"fCell",
    #"layerThickness0",
    #"maxLevelCell",
    #"minLevelCell",
    #"normalVelocity0",
    #"refBottomDepth",
    #"refZMid",
    #"restingThickness",
    #"salinity0",
    #"temperature0",
    #"vertCoordMovementWeights",
]
taglist_d = ["bottomDepth"]  # ["bed_elevation", "bottomDepth", "bottomDepthObserved", "fCell"]
taglist_i = ["maxLevelCell", "minLevelCell"]
taglist_z = ["refBottomDepth"]
taglist_3d = ["timeDaily_avg_ssh"]
#taglist_4d = ["temperature", "salinity", "layerThickness"]
taglist_4d = ["timeDaily_avg_activeTracers_temperature", "timeDaily_avg_activeTracers_salinity", "timeDaily_avg_layerThickness", "timeDaily_avg_velocityZonal", "timeDaily_avg_velocityMeridional" ]
#taglist_42d = ["temperature", "salinity"]
taglist_42d = []
hdf5_filename = "mpas_grid.h5m"
nc_filename = "data/dailyMeanOutput.0003-06-20.nc"
nc_gridfile = "data/init_data.nc"

taglist_i=[]
taglist_z=[]
taglist_4d=[]
taglist_42d = []

timestep=0
datathreshold = 1e30

opts, args = getopt.getopt(sys.argv[1:],"ht:m:d:g:",["timestep=", "moab=", "data=", "grid="])
for opt, arg in opts:
      if opt == '-h':
         print ('python convert_mpas_tags.py -d <MPAS Data file> -g <MPAS Grid file> -m <MOAB Outputfile> -t <timestep>')
         sys.exit()
      elif opt in ("-t", "--timestep"):
          timestep = int(arg)
      elif opt in ("-m", "--moab"):
         hdf5_filename = arg
      elif opt in ("-d", "--data"):
         nc_filename = arg
      elif opt in ("-g", "--grid"):
         nc_gridfile = arg

print ('Input MPAS grid file is  ', nc_gridfile)
print ('Input MPAS data file is  ', nc_filename)
print ('Output MPAS h5m file is  ', hdf5_filename)
print ('Current data timestep is ', timestep)

# start a MOAB instance
mb = core.Core()

# load the file
mb.load_file("mpas_grid_raw.h5m")

# get the root set of the MOAB instance
root_set = mb.get_root_set()

gidTag = mb.tag_get_handle("GLOBAL_ID")

# query the root set for all triangles
polys = mb.get_entities_by_type(root_set, types.MBPOLYGON, recur=True)
print("Found " + str(polys.size()) + " polygons in this model.")

# similar query for vertices
verts = mb.get_entities_by_type(root_set, types.MBVERTEX, recur=True)
print("Found " + str(verts.size()) + " vertices in this model.\n")

gids = mb.tag_get_data(gidTag, polys)

gids = gids.reshape((polys.size())) - 1

for etag in taglist_clean:
    print("Deleting", etag)
    thandle = mb.tag_get_handle(etag)
    mb.tag_delete(thandle)

tdata = np.zeros((polys.size()))
idata = np.zeros((polys.size()), dtype=np.int32)

# get the tag data
ncf = netCDF4.Dataset(nc_filename, "r")
nVertLevels = ncf.dimensions["nVertLevels"].size

ncg = netCDF4.Dataset(nc_gridfile, 'r')

if True:
    for dtag in taglist_d:
        print("\nAnalyzing", dtag)

        ncvar = ncg.variables[dtag]
        print(ncvar)

        # get the actual data out of the variable
        # tdata[gids[:]] = ncvar[:]
        tdata[:] = ncvar[gids[:]]

        # print(tdata[:50])

        print("Setting", dtag, "tag data")
        thandle = mb.tag_get_handle(dtag, 1, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, True)
        mb.tag_set_data(thandle, polys, tdata)

    for itag in taglist_i:
        print("\nAnalyzing", itag)

        ncvar = ncg.variables[itag]
        print(ncvar)

        # get the actual data out of the variable
        idata[:] = ncvar[gids[:]]

        print("Setting", itag, "tag data")
        thandle = mb.tag_get_handle(itag, 1, types.MB_TYPE_INTEGER, types.MB_TAG_DENSE, True)
        mb.tag_set_data(thandle, polys, idata)

    zdata = np.zeros((nVertLevels))
    for ztag in taglist_z:
        print("\nAnalyzing", ztag)

        ncvar = ncg.variables[ztag]
        print(ncvar)

        # get the actual data out of the variable
        zdata[:] = ncvar[:]

        print("Setting", ztag, "tag data")
        thandle = mb.tag_get_handle(
            ztag, nVertLevels, types.MB_TYPE_DOUBLE, types.MB_TAG_SPARSE, True
        )
        mb.tag_set_data(thandle, root_set, zdata)

    for fdtag in taglist_3d:
        print("\nAnalyzing", fdtag)

        ncvar = ncf.variables[fdtag]
        print(ncvar)

        # get the actual data out of the variable
        tdata[:] = ncvar[timestep, gids[:]]

        print("Setting", fdtag, "tag data")
        thandle = mb.tag_get_handle(fdtag,1,types.MB_TYPE_DOUBLE,types.MB_TAG_DENSE,True)
        mb.tag_set_data(thandle,polys,tdata)


    # tdata3d = np.zeros((polys.size(), nVertLevels))
    tdata3d = np.zeros((polys.size() * nVertLevels))

    for fdtag in taglist_4d:
        print("\nAnalyzing", fdtag)

        ncvar = ncf.variables[fdtag]
        print(ncvar)

        # get the actual data out of the variable
        tdata3d[:] = ncvar[timestep, gids[:], :].flatten()

        # reset data above threshold
        tdata3d[tdata3d > datathreshold] = 0

        # print(tdata3d[300*60:301*60])
        # print(tdata3d[301*60:302*60])

        print("Setting", fdtag, "tag data")
        thandle = mb.tag_get_handle(
            fdtag + "_3d", nVertLevels, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, True
        )
        mb.tag_set_data(thandle, polys, tdata3d)

    for fdtag in taglist_42d:
        print("\nAnalyzing", fdtag)

        ncvar = ncf.variables[fdtag]
        print(ncvar)

        # get the actual data out of the variable
        tdata[:] = ncvar[timestep, gids[:], 0].flatten()

        # reset data above threshold
        tdata[tdata > datathreshold] = 0

        print("Setting", fdtag, "tag data")
        thandle = mb.tag_get_handle(fdtag, 1, types.MB_TYPE_DOUBLE, types.MB_TAG_DENSE, True)
        mb.tag_set_data(thandle, polys, tdata)


## verify
if False:
    ncvar = ncf.variables["layerThickness"]
    print(ncvar)

    # get the actual data out of the variable
    tdata3d_a = ncvar[0, gids[1], :]
    tdata3d_b = ncvar[0, gids[2], :]
    tdata3d_c = ncvar[0, gids[3], :]

    print(tdata3d_a[:] - tdata3d_b[:])
    print(tdata3d_b[:] - tdata3d_c[:])

ncf.close()

print("Writing datasets in MOAB format to ", hdf5_filename) 
mb.write_file(hdf5_filename)

