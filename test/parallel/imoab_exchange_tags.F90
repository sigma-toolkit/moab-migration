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

         integer ierr, my_id, num_procs, pid, i, ix, iy, numv, nume
         integer dime, lco, mbtype, blockid, npe
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

         ! readopts( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS" )
         ! readoptsLnd( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION" )
         inputFileName = &
            MOAB_MESH_DIR &
            //'unittest/atm_c2x.h5m'//C_NULL_CHAR
         readopts = 'PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS'//C_NULL_CHAR
         num_layers = 2 ! two layer only

         call MPI_Init(ierr)
         call errorout(ierr, 'fail to initialize MPI')

         ierr = iMOAB_Initialize()
         call errorout(ierr, 'fail to initialize iMOAB')

         call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
         call errorout(ierr, 'fail to get MPI rank')

         call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
         call errorout(ierr, 'fail to get MPI size')

         ! find out MY process ID, and how many processes were started.
         if (my_id .eq. 0) then
            print *, " I'm process ", my_id, " out of ", &
               num_procs, " processes."
         end if

         !  give a component id
         compid = 7
         appname = 'IMTEST'//CHAR(0)
         ierr = iMOAB_RegisterApplication(appname, MPI_COMM_WORLD, &
                                          compid, pid)
         call errorout(ierr, 'fail to initialize fortran app')

         ierr = iMOAB_LoadMesh(pid, inputFileName, readopts, num_layers)
         call errorout(ierr, 'fail to load atm')

         !  see ghost elements
         dimgh = 2 ! will ghost quads, topological dim 2
         bridge_dim = 0 ! use vertex as bridge
         if (my_id .eq. 0) then
            print *, " Generating ", num_layers, " ghost layers"
         end if
         ierr = iMOAB_DetermineGhostEntities(pid, dimgh, num_layers, &
                                             bridge_dim)
         call errorout(ierr, 'fail to determine ghosts')

         tagTypes(:) = 1 ! DENSE_DOUBLE
         ierr = iMOAB_DefineTagStorage(pid, 'Sa_dens:Sa_pbot'//CHAR(0), tagTypes(1), 1, tagIndex(1))
         call errorout(ierr, 'failed to define the field tag Sa_dens ')

         entity_type(:) = 1 ! data is on elements (vertices = 0, elements = 1)
         ! let us now access the density field: Sa_dens
         ! and exchange tags to ensure the data is consistent
         ierr = iMOAB_SynchronizeTags(pid, 2, tagIndex, entity_type)

         ! write out the mesh file to disk, in parallel, if h5m
#ifdef MOAB_HAVE_HDF5
         write (outfile, '(A17,i1,A4)') 'ghost_exchanged_p', my_id, '.h5m'//CHAR(0)
         !outfile = 'ghost_exchanged.h5m'//CHAR(0)
         !wopts   = 'PARALLEL=WRITE_PART'//CHAR(0)
         wopts = CHAR(0)
         ierr = iMOAB_WriteMesh(pid, trim(outfile), wopts)
         call errorout(ierr, 'fail to write the mesh file')
#endif

         !  all done. de-register and finalize
         ierr = iMOAB_DeregisterApplication(pid)
         call errorout(ierr, 'fail to deregister application')

         ierr = iMOAB_Finalize()
         call errorout(ierr, 'fail to finalize iMOAB')

         call MPI_FINALIZE(ierr)
         call errorout(ierr, 'fail to finalize MPI')

      end

