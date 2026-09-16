!> @example DirectAccessF90.F90
!! \brief Direct (zero-copy) access to MOAB coordinate and tag storage from Fortran \n
!!
!! Every value a simulation exchanges with MOAB can be moved in two ways: copied through the
!! iMOAB accessors, or read and written in place through a pointer to MOAB's own storage.  This
!! example does the same computation both ways and checks that they agree bit for bit, so the
!! direct path is demonstrated and validated at the same time.
!!
!! The mesh is a 1d row of quads:
!!
!!  ----------------------
!!  |      |      |      |
!!  |      |      |      | ...
!!  |      |      |      |
!!  ----------------------
!!
!!    -#  Create the mesh through iMOAB_CreateVertices / iMOAB_CreateElements
!!    -#  Define tag1 and tag2 (3 doubles per quad) and vpe_tag (1 integer per quad)
!!    -#  Compute tag1, the quad centroids, using only the copying accessors
!!    -#  Compute tag2, the same centroids, writing through a direct tag pointer
!!    -#  Check tag1 == tag2 exactly
!!    -#  Translate every vertex in place through a direct coordinate pointer, then read the
!!        coordinates back with the copying accessor and check the shift landed exactly
!!    -#  Fill vpe_tag through a direct integer tag pointer and read it back with the copying
!!        accessor
!!
!! This supersedes DirectAccessNoHolesF90.F90, which reached MOAB's storage through iMesh.  The
!! "no holes" restriction is gone: the direct-access entry points report how many entities each
!! contiguous run covers, and the loops below simply continue until the whole range is consumed,
!! so a mesh whose handles are not one single block works unchanged.
!!
!! <b>To compile</b>: \n
!!    make DirectAccessF90 \n
!! <b>To run</b>: ./DirectAccessF90 \n
!!
program DirectAccessF90

   use iso_c_binding
   use iMOAB
   implicit none

#include "moab/MOABConfig.h"
#ifdef MOAB_HAVE_MPI
#  include "mpif.h"
#endif

   integer, parameter :: NQUADS = 1000
   integer, parameter :: NVERTS = 2*( NQUADS + 1 )
   integer, parameter :: VPE    = 4            ! vertices per quad
   integer, parameter :: MBQUAD = 3            ! moab::MBQUAD
   integer, parameter :: ON_ELEMS = 1   ! iMOAB entity_type; 0 would be the vertices
   double precision, parameter :: DX = 0.25d0, DY = -0.5d0, DZ = 2.0d0

   integer :: ierr, pid, compid, i, j, e, v, base, ndisagree
   integer :: nghlay, tagtype, ncomp, tagindex, blockid
   integer :: startidx, count, comps, nruns
   integer :: nvinfo(3), neinfo(3), nbinfo(3), nsbc(3), nvbc(3)
   integer :: vpe_in_block, nelem_in_block

   double precision :: coords(3*NVERTS), coords_after(3*NVERTS)
   double precision :: tag1(3*NQUADS), readback(3*NQUADS)
   integer :: conn(VPE*NQUADS), vpe_readback(NQUADS), blockids(16)

   ! Pointers into MOAB's own storage.  No C_PTR and no C_F_POINTER appear anywhere in this
   ! program: the iMOAB Fortran module hands back Fortran pointers directly.
   double precision, pointer :: xp(:), yp(:), zp(:)
   double precision, pointer :: tag2p(:)
   integer, pointer :: vpep(:)

#ifdef MOAB_HAVE_MPI
   call MPI_INIT( ierr )
#endif

   ierr = iMOAB_Initialize()
   call check( ierr, 'iMOAB_Initialize' )

   compid = 1
#ifdef MOAB_HAVE_MPI
   ierr = iMOAB_RegisterApplication( 'DIRECTACCESS', MPI_COMM_WORLD, compid, pid )
#else
   ierr = iMOAB_RegisterApplication( 'DIRECTACCESS', compid, pid )
#endif
   call check( ierr, 'iMOAB_RegisterApplication' )

   ! ---------------------------------------------------------------- build the mesh
   ! Vertices are numbered in layers: 2*i is the bottom of column i, 2*i+1 the top.
   do i = 0, NQUADS
      coords( 3*( 2*i ) + 1 )     = dble( i )
      coords( 3*( 2*i ) + 2 )     = 0.0d0
      coords( 3*( 2*i ) + 3 )     = 0.0d0
      coords( 3*( 2*i + 1 ) + 1 ) = dble( i )
      coords( 3*( 2*i + 1 ) + 2 ) = 1.0d0
      coords( 3*( 2*i + 1 ) + 3 ) = 0.0d0
   end do

   i = 3*NVERTS
   j = 3
   ierr = iMOAB_CreateVertices( pid, i, j, coords )
   call check( ierr, 'iMOAB_CreateVertices' )

   ! iMOAB_CreateElements takes 1-based vertex indices, unlike
   ! iMOAB_GetBlockElementConnectivities below, which returns 0-based ones.
   do i = 0, NQUADS - 1
      conn( VPE*i + 1 ) = 2*i + 1
      conn( VPE*i + 2 ) = 2*i + 3
      conn( VPE*i + 3 ) = 2*i + 4
      conn( VPE*i + 4 ) = 2*i + 2
   end do

   blockid = 100
   i = NQUADS
   j = MBQUAD
   e = VPE
   ierr = iMOAB_CreateElements( pid, i, j, e, conn, blockid )
   call check( ierr, 'iMOAB_CreateElements' )

   ierr = iMOAB_UpdateMeshInfo( pid )
   call check( ierr, 'iMOAB_UpdateMeshInfo' )

   ierr = iMOAB_GetMeshInfo( pid, nvinfo, neinfo, nbinfo, nsbc, nvbc )
   call check( ierr, 'iMOAB_GetMeshInfo' )
   if ( nvinfo(1) /= NVERTS .or. neinfo(1) /= NQUADS ) then
      print *, 'unexpected mesh size: ', nvinfo(1), ' vertices, ', neinfo(1), ' elements'
      call exit( 1 )
   end if

   ! iMOAB_GetBlockID is what builds the block-id lookup table; the calls below that take a
   ! block id fail with MB_FAILURE until it has run at least once.
   i = nbinfo(1)
   ierr = iMOAB_GetBlockID( pid, i, blockids )
   call check( ierr, 'iMOAB_GetBlockID' )

   ierr = iMOAB_GetBlockInfo( pid, blockid, vpe_in_block, nelem_in_block )
   call check( ierr, 'iMOAB_GetBlockInfo' )

   ! ---------------------------------------------------------------- define the tags
   tagtype = 1   ! dense, double
   ncomp = 3
   ierr = iMOAB_DefineTagStorage( pid, 'tag1', tagtype, ncomp, tagindex )
   call check( ierr, 'iMOAB_DefineTagStorage(tag1)' )
   ierr = iMOAB_DefineTagStorage( pid, 'tag2', tagtype, ncomp, tagindex )
   call check( ierr, 'iMOAB_DefineTagStorage(tag2)' )

   tagtype = 0   ! dense, integer
   j = 1
   ierr = iMOAB_DefineTagStorage( pid, 'vpe_tag', tagtype, j, tagindex )
   call check( ierr, 'iMOAB_DefineTagStorage(vpe_tag)' )

   ! ------------------------------------------------- tag1: the copying path
   i = 3*NVERTS
   ierr = iMOAB_GetVisibleVerticesCoordinates( pid, i, coords )
   call check( ierr, 'iMOAB_GetVisibleVerticesCoordinates' )

   i = VPE*NQUADS
   ierr = iMOAB_GetBlockElementConnectivities( pid, blockid, i, conn )
   call check( ierr, 'iMOAB_GetBlockElementConnectivities' )

   do e = 0, NQUADS - 1
      tag1( 3*e + 1 ) = 0.0d0
      tag1( 3*e + 2 ) = 0.0d0
      tag1( 3*e + 3 ) = 0.0d0
      do j = 0, VPE - 1
         v = conn( VPE*e + j + 1 )          ! 0-based local vertex index
         tag1( 3*e + 1 ) = tag1( 3*e + 1 ) + coords( 3*v + 1 )
         tag1( 3*e + 2 ) = tag1( 3*e + 2 ) + coords( 3*v + 2 )
         tag1( 3*e + 3 ) = tag1( 3*e + 3 ) + coords( 3*v + 3 )
      end do
      tag1( 3*e + 1 ) = tag1( 3*e + 1 ) / dble( VPE )
      tag1( 3*e + 2 ) = tag1( 3*e + 2 ) / dble( VPE )
      tag1( 3*e + 3 ) = tag1( 3*e + 3 ) / dble( VPE )
   end do

   i = 3*NQUADS
   ierr = iMOAB_SetDoubleTagStorage( pid, 'tag1', i, ON_ELEMS, tag1 )
   call check( ierr, 'iMOAB_SetDoubleTagStorage(tag1)' )

   ! ------------------------------------------------- tag2: the direct path
   ! The same arithmetic in the same order, but written straight into MOAB's storage rather
   ! than into a local array that iMOAB then copies.  start_index is 0-based, and count says
   ! how many elements the pointer covers, so the loop walks one contiguous run at a time.
   startidx = 0
   nruns = 0
   do while ( startidx < NQUADS )
      ierr = iMOAB_GetDoubleTagStoragePointer( pid, 'tag2', ON_ELEMS, startidx, count, comps, tag2p )
      call check( ierr, 'iMOAB_GetDoubleTagStoragePointer(tag2)' )
      call check_run( count, 'tag2' )
      nruns = nruns + 1
      if ( comps /= 3 ) then
         print *, 'tag2 has ', comps, ' components, expected 3'
         call exit( 1 )
      end if
      do i = 0, count - 1
         e = startidx + i                   ! element index within the whole local range
         base = comps*i                     ! offset within this run's pointer
         tag2p( base + 1 ) = 0.0d0
         tag2p( base + 2 ) = 0.0d0
         tag2p( base + 3 ) = 0.0d0
         do j = 0, VPE - 1
            v = conn( VPE*e + j + 1 )
            tag2p( base + 1 ) = tag2p( base + 1 ) + coords( 3*v + 1 )
            tag2p( base + 2 ) = tag2p( base + 2 ) + coords( 3*v + 2 )
            tag2p( base + 3 ) = tag2p( base + 3 ) + coords( 3*v + 3 )
         end do
         tag2p( base + 1 ) = tag2p( base + 1 ) / dble( VPE )
         tag2p( base + 2 ) = tag2p( base + 2 ) / dble( VPE )
         tag2p( base + 3 ) = tag2p( base + 3 ) / dble( VPE )
      end do
      startidx = startidx + count
   end do

   ! Read tag2 back the copying way, so the comparison exercises both directions.
   i = 3*NQUADS
   ierr = iMOAB_GetDoubleTagStorage( pid, 'tag2', i, ON_ELEMS, readback )
   call check( ierr, 'iMOAB_GetDoubleTagStorage(tag2)' )

   ndisagree = 0
   do i = 1, 3*NQUADS
      if ( differs( tag1( i ), readback( i ) ) ) ndisagree = ndisagree + 1
   end do
   if ( ndisagree /= 0 ) then
      print *, 'FAILED: tag1 and tag2 disagree in ', ndisagree, ' of ', 3*NQUADS, ' values'
      call exit( 1 )
   end if
   print *, 'tag written through a direct pointer matches the copied tag exactly, in', nruns, 'run(s)'

   ! ------------------------------------------------- direct coordinate access
   ! A streaming, in-place update: exactly the shape of a per-timestep field update, and the
   ! case where avoiding the copy matters most.
   startidx = 0
   nruns = 0
   do while ( startidx < NVERTS )
      ierr = iMOAB_GetVertexCoordinatesPointer( pid, startidx, count, xp, yp, zp )
      call check( ierr, 'iMOAB_GetVertexCoordinatesPointer' )
      call check_run( count, 'coordinates' )
      nruns = nruns + 1
      do i = 1, count
         xp( i ) = xp( i ) + DX
         yp( i ) = yp( i ) + DY
         zp( i ) = zp( i ) + DZ
      end do
      startidx = startidx + count
   end do

   i = 3*NVERTS
   ierr = iMOAB_GetVisibleVerticesCoordinates( pid, i, coords_after )
   call check( ierr, 'iMOAB_GetVisibleVerticesCoordinates after shift' )

   ndisagree = 0
   do v = 0, NVERTS - 1
      if ( differs( coords_after( 3*v + 1 ), coords( 3*v + 1 ) + DX ) ) ndisagree = ndisagree + 1
      if ( differs( coords_after( 3*v + 2 ), coords( 3*v + 2 ) + DY ) ) ndisagree = ndisagree + 1
      if ( differs( coords_after( 3*v + 3 ), coords( 3*v + 3 ) + DZ ) ) ndisagree = ndisagree + 1
   end do
   if ( ndisagree /= 0 ) then
      print *, 'FAILED: ', ndisagree, ' of ', 3*NVERTS, ' coordinates did not shift as expected'
      call exit( 1 )
   end if
   print *, 'coordinates written through a direct pointer read back exactly, in', nruns, 'run(s)'

   ! ------------------------------------------------- direct integer tag access
   startidx = 0
   nruns = 0
   do while ( startidx < NQUADS )
      ierr = iMOAB_GetIntTagStoragePointer( pid, 'vpe_tag', ON_ELEMS, startidx, count, comps, vpep )
      call check( ierr, 'iMOAB_GetIntTagStoragePointer(vpe_tag)' )
      call check_run( count, 'vpe_tag' )
      nruns = nruns + 1
      do i = 1, count
         vpep( i ) = vpe_in_block
      end do
      startidx = startidx + count
   end do

   i = NQUADS
   ierr = iMOAB_GetIntTagStorage( pid, 'vpe_tag', i, ON_ELEMS, vpe_readback )
   call check( ierr, 'iMOAB_GetIntTagStorage(vpe_tag)' )

   ndisagree = 0
   do i = 1, NQUADS
      if ( vpe_readback( i ) /= VPE ) ndisagree = ndisagree + 1
   end do
   if ( ndisagree /= 0 ) then
      print *, 'FAILED: ', ndisagree, ' of ', NQUADS, ' integer tag values are wrong'
      call exit( 1 )
   end if
   print *, 'integer tag written through a direct pointer reads back exactly, in', nruns, 'run(s)'

   ! ------------------------------------------------- a colon-separated list must be refused
   ! Each MOAB tag has its own allocation, so one pointer cannot describe two tags.  The API
   ! says so rather than quietly returning only the first.
   ierr = iMOAB_GetDoubleTagStoragePointer( pid, 'tag1:tag2', ON_ELEMS, 0, count, comps, tag2p )
   if ( ierr == 0 ) then
      print *, 'FAILED: a colon-separated tag list was accepted for direct access'
      call exit( 1 )
   end if
   print *, 'a colon-separated tag list is correctly rejected for direct access'

   ierr = iMOAB_DeregisterApplication( pid )
   call check( ierr, 'iMOAB_DeregisterApplication' )
   ierr = iMOAB_Finalize()
   call check( ierr, 'iMOAB_Finalize' )

#ifdef MOAB_HAVE_MPI
   call MPI_FINALIZE( ierr )
#endif

   print *, 'All direct-access results agree with the copying API, success!'

contains

   subroutine check( code, what )
      integer, intent(in) :: code
      character(len=*), intent(in) :: what
      if ( code /= 0 ) then
         print *, 'FAILED: ', what, ' returned ', code
         call exit( 1 )
      end if
   end subroutine check

   !> A run of zero entities would make the walking loops spin forever, so treat it as a failure
   !! rather than let the example hang.
   subroutine check_run( n, what )
      integer, intent(in) :: n
      character(len=*), intent(in) :: what
      if ( n <= 0 ) then
         print *, 'FAILED: direct access to ', what, ' returned a run of ', n, ' entities'
         call exit( 1 )
      end if
   end subroutine check_run

   !> Exact comparison is deliberate, and is the whole point of the checks: both paths run the
   !! same arithmetic in the same order over the same inputs, so any difference at all means
   !! the direct pointer is not aliasing the storage the copying accessor reads.
   logical function differs( a, b )
      double precision, intent(in) :: a, b
      differs = ( a /= b )
   end function differs

end program DirectAccessF90
