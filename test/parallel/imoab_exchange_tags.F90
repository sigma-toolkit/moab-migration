!
!      This program shows how to load a mesh in parallel into MOAB using iMOAB, with sufficient
!      information to resolve boundary sharing and exchange a layer of ghost information.
!
!       After resolving the sharing, we obtain a tag handle and exchange the data on ghost layers.
!       Then the individual mesh+tag data is saved to a file.
!
!       By default, this test is run on 2 processors
!
      SUBROUTINE errorout(ierr, message)
         integer ierr
         character*(*) message
         if (ierr .ne. 0) then
            print *, message
            call exit(1)
         end if
         return
      end

#include "moab/MOABConfig.h"
      program imoabexchangepartag_test

         use iMOAB
         implicit none

#include "mpif.h"

         integer ierr, rank, nprocs, pid
         integer compid
         character :: appname*10
         character(:), allocatable :: inputFileName
         character(:), allocatable :: readopts
         character :: outfile*100, wopts*100
!      used for ghosting
         integer dimgh, bridge_dim, num_layers
         integer, dimension(2) ::  entity_type
         integer, dimension(2) ::  tagIndex
         integer, dimension(2) :: tagTypes!  { DENSE_DOUBLE, DENSE_DOUBLE }
         integer nverts(3), nelem(3), nblocks(3), nsbc(3), ndbc(3)

         ! readopts( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS" )
         ! readoptsLnd( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION" )
         !inputFileName = &
         !   MOAB_MESH_DIR &
         !   //'unittest/atm_c2x.h5m'//C_NULL_CHAR
         !readopts = 'PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;' &
         !           //'PARALLEL_RESOLVE_SHARED_ENTS'//C_NULL_CHAR
         inputFileName = MOAB_MESH_DIR &
              // 'unittest/io/domain.ocn.ne4np4_oQU240.160614.nc'//C_NULL_CHAR
         readopts = 'PARALLEL=READ_PART;PARTITION_METHOD=SQIJ;VARIABLE=;' &
                    //'REPARTITION'//C_NULL_CHAR

         num_layers = 0 ! initialize to no ghost layer for now

         call MPI_Init(ierr)
         call errorout(ierr, 'fail to initialize MPI')

         ierr = iMOAB_Initialize()
         call errorout(ierr, 'fail to initialize iMOAB')

         call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
         call errorout(ierr, 'fail to get MPI rank')

         call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
         call errorout(ierr, 'fail to get MPI size')

         ! find out MY process ID, and how many processes were started.
         if (rank .eq. 0) then
            print *, " I'm process ", rank, " out of ", &
               nprocs, " processes."
         end if

         !  give a component id
         compid = 100
         appname = 'HaloExchangeApp'//CHAR(0)

         ! first, let us register the application
         ierr = iMOAB_RegisterApplication(appname, MPI_COMM_WORLD, &
                                          compid, pid)
         call errorout(ierr, 'fail to initialize fortran app')

         ! load the mesh in parallel with specified options
         ierr = iMOAB_LoadMesh(pid, inputFileName, readopts, num_layers)
         call errorout(ierr, 'fail to load atm')

         !  see ghost elements
         dimgh = 2 ! will ghost quads, topological dim 2
         num_layers = 2 ! let us get two ghost layers
         bridge_dim = 0 ! use vertex as bridge
         if (rank .eq. 0) then
            print *, " Generating ", num_layers, " ghost layers"
         end if
         ierr = iMOAB_DetermineGhostEntities(pid, dimgh, num_layers, &
                                             bridge_dim)
         call errorout(ierr, 'fail to determine ghosts')

         ! let us get some information about the partitioned mesh and print
         ierr = iMOAB_GetMeshInfo(pid, nverts, nelem, nblocks, nsbc, ndbc)
         call errorout(ierr, 'Error: failed to get mesh info ')
         if (rank .eq. 0) then
            print *,  "MOAB vertices: owned=", nverts(1), &
                              ", ghosted=", nverts(2), ", total=", nverts(3)
            print *,  "MOAB elements: owned=", nelem(1), &
                              ", ghosted=", nelem(2), ", total=", nelem(3)
         endif

         tagTypes(:) = 1 ! DENSE_DOUBLE
         ierr = iMOAB_DefineTagStorage(pid, 'Sa_dens:Sa_pbot'//CHAR(0), &
                                       tagTypes(1), 1, tagIndex(1))
         call errorout(ierr, 'failed to define the field tag Sa_dens ')

         entity_type(:) = 1 ! data is on elements (vertices = 0, elements = 1)
         ! let us now synchronize the density/pressure fields: Sa_dens, Sa_pbot
         ! and exchange tag values to ensure the data is consistent in the ghost
         ! halo regions
         ierr = iMOAB_SynchronizeTags(pid, 2, tagIndex, entity_type)

         ! write out the mesh file to disk, in parallel, if h5m
#ifdef MOAB_HAVE_HDF5
         outfile = 'ghost_exchanged'//CHAR(0)
         ierr = iMOAB_WriteLocalMesh(pid, trim(outfile))
         call errorout(ierr, 'fail to write the local mesh files')
#endif

         !  all done. de-register and finalize
         ierr = iMOAB_DeregisterApplication(pid)
         call errorout(ierr, 'fail to deregister application')

         ierr = iMOAB_Finalize()
         call errorout(ierr, 'fail to finalize iMOAB')

         call MPI_FINALIZE(ierr)
         call errorout(ierr, 'fail to finalize MPI')

      end

