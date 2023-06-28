from pymoab import core
from pymoab import types
from driver import test_driver_parallel, CHECK_EQ, CHECK_NOT_EQ, CHECK_ITER_EQ

from mpi4py import MPI
from pymoab import parallelcomm

import numpy as np
import os

bytes_per_char_ = np.array(["a"]).nbytes

def test_parallel_rank_size(mpicomm):

    mpirank = mpicomm.Get_rank()
    mpisize = mpicomm.Get_size()

    mb = core.Core()
    pc = parallelcomm.ParallelComm(mb, comm=mpicomm)

    pid = pc.get_id()
    rank = pc.rank()
    size = pc.size()

    assert pid >= 0
    CHECK_EQ(mpirank, rank)
    CHECK_EQ(mpisize, size)

def test_parallel_load_mesh(mpicomm):
    mb = core.Core()
    try:
        mb.load_file("parallel_file.h5m", file_set = None, readopts = 'PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS')
    except:
        try:
            print( """
            WARNING: .h5m file load failed. If hdf5 support is enabled in this
            build there could be a problem.
            """)
            mb.load_file("parallel_file.vtk", file_set = 0, readopts = '')
        except:
            raise(IOError, "Failed to load MOAB file.")

    #load into file_set
    mb1 = core.Core()
    file_set = mb1.create_meshset()
    try:
        mb1.load_file("parallel_file.h5m", file_set, readopts = 'PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS')
    except:
        try:
            print("""
            WARNING: .h5m parallel file load failed. If HDF5 support is enabled in this
            build there could be a problem.
            """)
            mb1.load_file("parallel_file.vtk",file_set, readopts = '')
        except:
            raise(IOError, "Failed to parallel load MOAB file.")

    ents = mb1.get_entities_by_type(file_set,types.MBMAXTYPE)
    CHECK_NOT_EQ(len(ents),0)

def test_parallel_write_mesh(mpicomm):
    mb = core.Core()
    mb.create_vertices(np.ones(3))

    try:
        mb.write_file("poutfile.h5m", output_sets = None, writeopts = 'PARALLEL=WRITE_PART')
        assert os.path.isfile("poutfile.h5m")
    except:
        try:
            print("""
            WARNING: .h5m file write failed. If HDF5 support is enabled in this
            build there could be a problem.
            """)
            mb.write_file("poutfile.vtk", output_sets = None, writeopts = '')
            assert os.path.isfile("outfile.vtk")
        except:
            raise(IOError, "Failed to parallel write MOAB file.")


def test_parallel_write_tags(mpicomm):
    """
    Test write tag functionality
    """

    # test values
    outfile = "parallel_write_tag_test.h5m"

    mb = core.Core()
    pc = parallelcomm.ParallelComm(mb, comm=mpicomm)
    rank = pc.rank()

    ##                  (1,1)3
    ##
    ##      PART 0               PART 1
    ##
    ## (0,0)1           (1,0)2              (2,0)4
    if rank == 0:
        coords = np.array((0,0,0,1,0,0,1,1,0),dtype='float64')
    else:
        coords = np.array((1,0,0,2,0,0,1,1,0),dtype='float64')
    vertices = mb.create_vertices(coords)
    CHECK_EQ(len(vertices),3)

    #create elements
    verts = np.array(((vertices[0],vertices[1],vertices[2]),),dtype='uint64')
    tris = mb.create_elements(types.MBTRI,verts)
    CHECK_EQ(len(tris),1)

    #check that the element is there via GLOBAL_ID tag
    global_id_tag = mb.tag_get_handle(types.GLOBAL_ID_TAG_NAME)
    if rank == 0:
        gids = np.array((1, 2, 3), dtype='int64')
    else:
        gids = np.array((2, 4, 3), dtype='int64')
    mb.tag_set_data(global_id_tag, vertices, gids)

    pset = pc.create_part()
    mb.add_entities(pset, tris)

    pc.resolve_shared_entities(pset, tris, resolve_dim=2, shared_dim=1)

    # pc.resolve_shared_entities(pset, 2, 1)

    # pc.assign_global_ids(setid=pset, dimension=2, startid=1, largestdimonly=False, isparallel=True, ownedonly=False)
    # rank = np.array((pc.rank()), dtype='int64')
    # mb.tag_set_data(pc.partition_tag(), pset, rank)

    # create writing tag
    write_tag = mb.tag_get_handle("WRITE",
                                  3,
                                  types.MB_TYPE_DOUBLE,
                                  types.MB_TAG_DENSE,
                                  create_if_missing=True)
    # set some data on that tag
    if rank == 0:
        data = [0.7071, 0.7071, 0.0, 0.5, 0.0, 0.5, -1.0, 1.0, 0.7071]
    else:
        data = [0.5, 0.0, 0.5, 0.7071, 0.7071, 0.0, -1.0, 1.0, 0.7071]

    mb.tag_set_data(write_tag, vertices, data)

    mb.write_file(outfile, output_tags = [write_tag], writeopts = 'PARALLEL=WRITE_PART')

    mb2 = core.Core()
    # mb2.load_file(outfile, readopts = 'PARALLEL=BCAST_DELETE;PARTITION=PARALLEL_TRIVIAL;PARALLEL_RESOLVE_SHARED_ENTS')
    mb2.load_file(outfile, readopts = 'PARALLEL=BCAST_DELETE;PARTITION=TRIVIAL;PARALLEL_RESOLVE_SHARED_ENTS')

    # pc.resolve_shared_entities(pset, tris, resolve_dim=2, shared_dim=1)

    vs = mb2.get_entities_by_type(0, types.MBVERTEX)

    CHECK_EQ(len(vs), 3)

    # get the write tag
    new_write_tag = mb2.tag_get_handle("WRITE")

    # make sure we can still get data for the write tag
    d = mb2.tag_get_data(new_write_tag, vs)

    # check if the data we set on vertices are correct
    if rank == 0:
        CHECK_ITER_EQ(d[0], [0.7071, 0.7071, 0.0])
        CHECK_ITER_EQ(d[1], [0.5, 0.0, 0.5])
        CHECK_ITER_EQ(d[2], [-1.0, 1.0, 0.7071])
    else:
        CHECK_ITER_EQ(d[0], [0.5, 0.0, 0.5])
        CHECK_ITER_EQ(d[1], [-1.0, 1.0, 0.7071])
        CHECK_ITER_EQ(d[2], [0.7071, 0.7071, 0.0])


def test_assign_global_ids(mpicomm):
    """
    Test assignment of global ID numbers
    """

    # test values
    outfile = "parallel_write_tag_test.h5m"

    mb = core.Core()
    pc = parallelcomm.ParallelComm(mb, comm=mpicomm)
    rank = pc.rank()

    ##                  (1,1)3
    ##
    ##      PART 0               PART 1
    ##
    ## (0,0)1           (1,0)2              (2,0)4
    if rank == 0:
        coords = np.array((0,0,0,1,0,0,1,1,0),dtype='float64')
    else:
        coords = np.array((1,0,0,2,0,0,1,1,0),dtype='float64')
    vertices = mb.create_vertices(coords)
    CHECK_EQ(len(vertices),3)

    #create elements
    verts = np.array(((vertices[0],vertices[1],vertices[2]),),dtype='uint64')
    tris = mb.create_elements(types.MBTRI,verts)
    CHECK_EQ(len(tris),1)

    #check that the element is there via GLOBAL_ID tag
    global_id_tag = mb.tag_get_handle(types.GLOBAL_ID_TAG_NAME)
    if rank == 0:
        gids = np.array((1, 2, 3), dtype='int64')
    else:
        gids = np.array((2, 4, 3), dtype='int64')
    mb.tag_set_data(global_id_tag, vertices, gids)

    egids = [[(rank+1)%2+1]] # rank 0: gid=2, rank 1: gid=1
    mb.tag_set_data(global_id_tag, tris, egids)

    pset = pc.create_part()
    mb.add_entities(pset, tris)

    pc.resolve_shared_entities(pset, tris, resolve_dim=2, shared_dim=1)

    pc.assign_global_ids(setid=pset, dimension=2, startid=1, largestdimonly=True, isparallel=True, ownedonly=False)

    e2gids = mb.tag_get_data(global_id_tag, tris)

    CHECK_EQ(e2gids, [[rank+1]])


if __name__ == "__main__":
    tests = [test_parallel_rank_size,
             test_parallel_load_mesh,
             test_parallel_write_mesh,
             test_parallel_write_tags,
             test_assign_global_ids
             ]

    test_driver_parallel(tests, MPI.COMM_WORLD)
