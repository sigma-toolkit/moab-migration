!This program shows how to do perform FV-FV projections in multiple ways in Fortran90
!
!The program workflows is as follows :
! - Setup a source( ATM ) and a target( OCN ) mesh in some PEs
! - Next migrate the meshes to the coupler PEs and compute source coverage mesh
! - Then compute the intersection mesh between the source and target meshes
! - Next compute the weights for the projection for 4 different methods
!    * First order projection
!    * Bilinear projection
!    * Second order projection
!    * Second order projection with CAAS operator
! - Then apply the weights to the source mesh and get the projected fields
!   for all the different methods
! - Finally send all the projected fields back to the target PEs
!
!   The program uses iMOAB Fortran90 interface to MOAB

SUBROUTINE errorout(ierr, message)
   integer ierr
   character*(*) message
   if (ierr .ne. 0) then
      print *, message
      call exit(1)
   end if
   return
end
!
SUBROUTINE check_baseline(baseline_file, nsize, gids, values, eps, rank, ierr)
   integer :: ierr
   character(*) ::  baseline_file
   integer :: nsize
   integer  :: gids (nsize)
   double precision  :: values(nsize)
   integer :: rank
   double precision :: eps 
   integer , allocatable :: allgids(:)
   double precision , allocatable :: allvals(:)
   integer :: unit, n ! n is for number of rows in the file
   
   
   unit = 21 + rank ! to differentiate them
   open(unit, file = baseline_file,status="old",action="read")
   ierr = 0
   n = 0
   do
      read(unit,*,end=1)
      n = n+1
   end do
1  rewind(unit)
   allocate(allgids(n))
   allocate(allvals(n))
   
   do i = 1,n
      read(unit,*) allgids(i), allvals(i)  ! we should have allgids from 1 to n, actually
   enddo
   
   do i = 1, nsize
      if ( abs( values(i) - allvals(gids(i)) ) .gt. eps) then
          print *, 'rank:', rank, ' index i', i, ' values:', values(i), &
           'gids:', gids(i),  'allvals(gids(i)): ',    allvals(gids(i))
          ierr = 1
      endif
   end do
    
   return
end

#include "moab/MOABConfig.h"

#define VERBOSE

#ifndef MOAB_MESH_DIR
#error Specify MOAB_MESH_DIR path
#endif

program imoab_coupler_fortran

   use iso_c_binding
   use iMOAB
   implicit none

#include "mpif.h"
#include "moab/MOABConfig.h"
   integer :: m ! for number of arguments ; if less than 1, exit
   integer :: global_comm
   integer :: rankInGlobalComm, numProcesses
   integer :: ierr
   integer :: my_id, num_procs
   integer :: jgroup   ! group for global comm
   character(:), allocatable :: atmFileName
   character(:), allocatable :: ocnFileName
   character(:), allocatable :: readopts, fileWriteOptions
   character(:), allocatable :: base_file3 !  baseline for bilinear example
   character(:), allocatable :: base_file4 !  baseline for higher order example
   character(:), allocatable :: base_file5 !  baseline for higher order example caas
   character :: appname*128
   character(30) :: nproc
   character(:), allocatable :: weights_identifier1
   character(:), allocatable :: disc_methods1, disc_methods2, dof_tag_names1, dof_tag_names2
   integer :: disc_orders1, disc_orders2
   character(:), allocatable :: atmocn_map_file_name, intx_from_file_identifier

   ! all groups and comms are on the same number of processes, for simplicity
   integer :: atmGroup, atmComm, ocnGroup, ocnComm, cplGroup, cplComm
   integer :: atmCouComm, ocnCouComm

   integer :: cmpatm, cplatm, cmpocn, cplocn, atmocnid ! global IDs
   integer :: cmpAtmPID, cplAtmPID ! iMOAB app ids
   integer :: cmpOcnPID, cplOcnPID ! iMOAB app ids
   integer :: cplAtmOcnPID ! intx pid
   integer :: nghlay, partScheme, context_id, nghlay_tgt
   integer :: fNoBubble, fMonotoneTypeID, fVolumetric, fNoConserve, fValidate, fInverseDistanceMap
   integer :: filter_type

   integer, dimension(2) ::  tagIndex
   integer, dimension (2) :: tagTypes!  { DENSE_DOUBLE, DENSE_DOUBLE }
   integer :: atmCompNDoFs ! = disc_orders[0] * disc_orders[0],
   integer :: ocnCompNDoFs !  = 1 /*FV*/
   character(:), allocatable :: fields, projectedFields, projectedFieldsBilin, projectedFieldsSecond, projectedFieldsCAAS
   integer, dimension(3) ::  nverts, nelem, nblocks, nsbc, ndbc
   double precision, allocatable :: vals(:) ! to set the double values to 0
   integer, allocatable :: gids(:) ! integers global ids
   integer :: i ! for loops
   integer :: storLeng, eetype ! for tags defs
   character(:), allocatable :: transferFields, outputFileOcn
   integer :: tagIndexIn2 ! not really needed
   integer :: type1, type2 ! for comm graph between atm cpl and atm coverage on cpl for ocean
   double precision  :: eps
   integer :: gnomonic

   cmpatm = 5
   cplatm = 6
   cmpocn = 17
   cplocn = 18
   atmocnid = 618

   call mpi_init(ierr)
   call errorout( ierr, 'mpi_init' )
   call mpi_comm_rank( MPI_COMM_WORLD, my_id, ierr )
   call errorout( ierr, 'fail to get MPI rank' )

   call mpi_comm_size( MPI_COMM_WORLD, num_procs, ierr )
   call errorout( ierr, 'fail to get MPI size' )
   call mpi_comm_dup( MPI_COMM_WORLD, global_comm, ierr )
   call errorout( ierr, 'fail to get global comm duplicate' )

   call MPI_Comm_group( global_comm, jgroup, ierr );  !  all processes in jgroup
   call errorout(ierr, 'fail to get joint group')
   ! readopts( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS" )
   ! readoptsLnd( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION" )
   atmFileName = &
     MOAB_MESH_DIR &
     //'unittest/atm_c2x.h5m'//C_NULL_CHAR
   ocnFileName = &
     MOAB_MESH_DIR &
     //'unittest/wholeOcn.h5m'//C_NULL_CHAR
     
   base_file3 = &
     MOAB_MESH_DIR &
     //'unittest/baseline3.txt'//C_NULL_CHAR
   base_file4 = &
     MOAB_MESH_DIR &
     //'unittest/baseline4.txt'//C_NULL_CHAR
   base_file5 = &
     MOAB_MESH_DIR &
     //'unittest/baseline5.txt'//C_NULL_CHAR
     

   ! all comms span the whole world, for simplicity
   atmComm = MPI_COMM_NULL
   ocnComm = MPI_COMM_NULL
   cplComm = MPI_COMM_NULL
   atmCouComm = MPI_COMM_NULL
   ocnCouComm = MPI_COMM_NULL
   call mpi_comm_dup(global_comm, atmComm, ierr)
   call mpi_comm_dup( global_comm, ocnComm, ierr )
   call mpi_comm_dup( global_comm, cplComm, ierr )
   call mpi_comm_dup( global_comm, atmCouComm, ierr )
   call mpi_comm_dup( global_comm, ocnCouComm, ierr )

   ! all groups of interest are easy breezy
   call MPI_Comm_group( atmComm, atmGroup, ierr )
   call MPI_Comm_group( ocnComm, ocnGroup, ierr )
   call MPI_Comm_group( cplComm, cplGroup, ierr )

   readopts ='PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS'//C_NULL_CHAR
   nghlay = 0 ! no ghost layers
#ifdef MOAB_HAVE_ZOLTAN
   partScheme = 2  ! RCB with zoltan
#else
   partScheme = 0  ! Trivial partitioner
#endif

   if (my_id .eq. 0) then
      print *, ' number of tasks: ', num_procs
      print *, ' Atm file: ', atmFileName
      print *, ' Ocn file: ', ocnFileName
      print *, ' using partitioner: ', partScheme
   end if

   ierr = iMOAB_Initialize()
   appname = 'ATM'//C_NULL_CHAR
   ierr = iMOAB_RegisterApplication(appname, atmComm, cmpatm, cmpatmPid)
   appname = 'ATMX'//C_NULL_CHAR
   ierr = iMOAB_RegisterApplication(appname, cplComm, cplatm, cplatmPid)
   appname = 'OCN'//C_NULL_CHAR
   ierr = iMOAB_RegisterApplication(appname, ocnComm, cmpocn, cmpocnPid)
   appname = 'OCNX'//C_NULL_CHAR
   ierr = iMOAB_RegisterApplication(appname, cplComm, cplocn, cplocnPid)

   appname = 'ATMOCN'//C_NULL_CHAR
   ierr = iMOAB_RegisterApplication(appname, cplComm, atmocnid, cplAtmOcnPID)

   ! read atm and migrate
   if (atmComm .NE. MPI_COMM_NULL) then
      ierr = iMOAB_LoadMesh(cmpatmPid, atmFileName, readopts, nghlay)
      call errorout(ierr, 'fail to load atm')
      ierr = iMOAB_SendMesh(cmpatmPid, atmCouComm, cplGroup, cplatm, partScheme)
      call errorout(ierr, 'fail to send atm')
   end if
   if (cplComm .NE. MPI_COMM_NULL) then
      ierr = iMOAB_ReceiveMesh(cplatmPid, atmCouComm, atmGroup, cmpatm) !
      call errorout(ierr, 'fail to receive atm')
   end if

   ! we can now free the sender buffers
   if (atmComm .NE. MPI_COMM_NULL) then
      context_id = cplatm
      ierr = iMOAB_FreeSenderBuffers(cmpatmPid, context_id)
      call errorout(ierr, 'fail to free atm buffers')
   end if

   ! read ocn and migrate
   if (ocnComm .NE. MPI_COMM_NULL) then
      ierr = iMOAB_LoadMesh(cmpocnPid, ocnFileName, readopts, nghlay)
      call errorout(ierr, 'fail to load ocn')
      ierr = iMOAB_SendMesh(cmpocnPid, ocnCouComm, cplGroup, cplocn, partScheme)
      call errorout(ierr, 'fail to send ocn')
   end if
   if (cplComm .NE. MPI_COMM_NULL) then
      ierr = iMOAB_ReceiveMesh(cplocnPid, ocnCouComm, ocnGroup, cmpocn) !
      call errorout(ierr, 'fail to receive ocn')
   end if

   ! we can now free the sender buffers
   if (ocnComm .NE. MPI_COMM_NULL) then
      context_id = cplocn
      ierr = iMOAB_FreeSenderBuffers(cmpocnPid, context_id)
      call errorout(ierr, 'fail to free ocn buffers')
   end if

   if (cplComm .NE. MPI_COMM_NULL) then

      ! set the ghost layers on the coupler for the source mesh
      nghlay = 1 ! number of ghost layers
      nghlay_tgt = 0
      ierr = iMOAB_SetMapGhostLayers( cplAtmOcnPID, nghlay, nghlay_tgt )
      call errorout(ierr, 'failed to set number of ghost layers on ATM mesh')

      ierr = iMOAB_ComputeMeshIntersectionOnSphere(cplAtmPID, cplOcnPID, cplAtmOcnPID)
      ! coverage mesh was computed here, for cplAtmPID, atm on coupler pes
      ! basically, atm was redistributed according to target (ocean) partition, to "cover" the
      !ocean partitions check if intx valid, write some h5m intx file
      call errorout(ierr, 'cannot compute intersection')

#ifdef VERBOSE
      ierr = iMOAB_WriteLocalMesh(cplAtmOcnPID, 'intx_ao')
      call errorout(ierr, 'could not write intersection mesh to disk')
#endif

   end if

   if (cplComm .NE. MPI_COMM_NULL) then
      ! the new graph will be for sending data from atm cpl  to atm coverage mesh for ocean
      ! it involves cpl atm app;  cplAtmPID
      ! results are in cplAtmOcnPID, intx mesh; remapper also has some info about coverage mesh
      ! after this, the sending of tags from atm pes to coupler pes will use the new par comm
      ! graph, that has more precise info about what to send for ocean cover ; every time, we
      ! will use the element global id, which should uniquely identify the element
      type1 = 3;
      type2 = 3;
      ierr = iMOAB_ComputeCommGraph( cplAtmPID, cplAtmOcnPID, cplComm, cplGroup, cplGroup, &
         type1, type2, cplatm, atmocnid )
      call errorout(ierr, 'cannot recompute direct coverage graph for atm coverage for ocean')
   end if

   weights_identifier1 = 'scalar'//C_NULL_CHAR
   disc_methods1 = 'fv'//C_NULL_CHAR
   disc_methods2 = 'fv'//C_NULL_CHAR
   disc_orders1 = 1
   disc_orders2 = 1
   dof_tag_names1 = 'GLOBAL_ID'//C_NULL_CHAR
   dof_tag_names2 = 'GLOBAL_ID'//C_NULL_CHAR
   ! fMonotoneTypeID = 0, fVolumetric = 0, fValidate = 1, fNoConserve = 0, fNoBubble = 1
   fNoBubble = 1
   fMonotoneTypeID = 0
   fVolumetric = 0
   fNoConserve = 0
   fValidate = 0
   fInverseDistanceMap = 0
   filter_type = 0

   if (cplComm .NE. MPI_COMM_NULL) then

      ierr = iMOAB_ComputeScalarProjectionWeights( &
             cplAtmOcnPID, weights_identifier1, disc_methods1, disc_orders1, &
             disc_methods2, disc_orders2, ""//C_NULL_CHAR, fNoBubble, fMonotoneTypeID, fVolumetric, &
             fInverseDistanceMap, fNoConserve, &
             fValidate, dof_tag_names1, dof_tag_names2)
      call errorout(ierr, 'cannot compute scalar first order projection weights')

      ierr = iMOAB_ComputeScalarProjectionWeights( &
             cplAtmOcnPID, "bilinear"//C_NULL_CHAR, disc_methods1, 1, &
             disc_methods2, 1, "bilin"//C_NULL_CHAR, fNoBubble, fMonotoneTypeID, fVolumetric, &
             fInverseDistanceMap, fNoConserve, &
             fValidate, "GLOBAL_ID"//C_NULL_CHAR, "GLOBAL_ID"//C_NULL_CHAR)
      call errorout(ierr, 'cannot compute scalar bilinear projection weights')

      ierr = iMOAB_ComputeScalarProjectionWeights( &
             cplAtmOcnPID, "secondorder"//C_NULL_CHAR, disc_methods1, 2, &
             disc_methods2, 2, ""//C_NULL_CHAR, fNoBubble, fMonotoneTypeID, fVolumetric, &
             fInverseDistanceMap, fNoConserve, &
             fValidate, "GLOBAL_ID"//C_NULL_CHAR, "GLOBAL_ID"//C_NULL_CHAR)
      call errorout(ierr, 'cannot compute scalar 2nd order projection weights')

#if defined( MOAB_HAVE_NETCDF ) && defined( VERBOSE )
      write(nproc,"(I0.2)")num_procs !
      atmocn_map_file_name = 'atm_ocn_map_second_n'//trim(nproc)//'.nc'//C_NULL_CHAR
      ierr = iMOAB_WriteMappingWeightsToFile( cplAtmOcnPID, "secondorder"//C_NULL_CHAR, atmocn_map_file_name)
      call errorout(ierr, 'failed to write map file to disk')
#endif

   end if

   ! start copy
   tagTypes(1) = 1 ! somehow, DENSE_DOUBLE give 0, while it should be 1; maybe moab::DENSE_DOUBLE ?
   tagTypes(2) = 1 ! ! DENSE_DOUBLE
   atmCompNDoFs = disc_orders1*disc_orders1
   ocnCompNDoFs = 1 ! /*FV*/

   fields = 'Sa_dens:Sa_pbot'//C_NULL_CHAR
   projectedFields = 'Sa_dens_proj:Sa_pbot_proj'//C_NULL_CHAR
   projectedFieldsBilin = 'Sa_dens_bilin_proj:Sa_pbot_bilin_proj'//C_NULL_CHAR
   projectedFieldsSecond = 'Sa_dens_o2_proj:Sa_pbot_o2_proj'//C_NULL_CHAR
   projectedFieldsCAAS = 'Sa_dens_o2_caas_proj:Sa_pbot_o2_caas_proj'//C_NULL_CHAR
   transferFields = 'Sa_dens_proj:Sa_pbot_proj:Sa_dens_bilin_proj:Sa_pbot_bilin_proj:Sa_dens_o2_proj:'//&
                    'Sa_pbot_o2_proj:Sa_dens_o2_caas_proj:Sa_pbot_o2_caas_proj'//C_NULL_CHAR

   if (cplComm .NE. MPI_COMM_NULL) then
      ierr = iMOAB_DefineTagStorage(cplAtmPID, fields, tagTypes(1), atmCompNDoFs, tagIndex(1))
      call errorout(ierr, 'failed to define the field tags a2oTbot:a2oUbot:a2oVbot ')
      ierr = iMOAB_DefineTagStorage(cplOcnPID, transferFields, tagTypes(2), ocnCompNDoFs, tagIndex(2))
      call errorout(ierr, 'failed to define the field tags a2oTbot_proj:a2oUbot_proj:a2oVbot_proj')
   end if

   if (atmComm .NE. MPI_COMM_NULL) then
      ! As always, use nonblocking sends
      ! this is for projection to ocean:
      ierr = iMOAB_SendElementTag(cmpAtmPID, fields, atmCouComm, cplatm)
      call errorout(ierr, 'cannot send tag values')

   end if

   if (cplComm .NE. MPI_COMM_NULL) then
      !// receive on atm on coupler pes, that was redistributed according to coverage
      ierr = iMOAB_ReceiveElementTag(cplAtmPID, fields, atmCouComm, cmpatm)
      call errorout(ierr, 'cannot receive tag values')
   end if

   ! we can now free the sender buffers
   if (atmComm .NE. MPI_COMM_NULL) then
      ierr = iMOAB_FreeSenderBuffers(cmpAtmPID, cplatm) !context is for ocean
      call errorout(ierr, 'cannot free buffers used to resend atm tag towards the coverage mesh')
   end if

   if (cplComm .ne. MPI_COMM_NULL) then

      outputFileOcn = "AtmOnCplF.h5m"//C_NULL_CHAR
      fileWriteOptions = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
      ierr = iMOAB_WriteMesh(cplAtmPID, outputFileOcn, fileWriteOptions)
      call errorout(ierr, 'could not write AtmOnCpl.h5m to disk')

   end if
   
    ! we need a second hop, to send from cpl atm to atm coverage for ocn 
    ! second hop, is from atm towards ocean, on coupler
    !  it should send from each part on coupler towards the coverage set that should form the
    ! rings around target cells (ocean)
    ! basically we should send to more cells than needed just for intersection
    !  TODO
    if( cplComm .ne. MPI_COMM_NULL ) then
        ! send using the par comm graph computed by iMOAB_ComputeCommGraph
        ierr = iMOAB_SendElementTag( cplAtmPID, fields, cplComm, atmocnid )
        call errorout( ierr, "cannot send tag values towards coverage mesh for bilinear map" )

        ierr = iMOAB_ReceiveElementTag( cplAtmOcnPID, fields, cplComm, cplatm )
        call errorout( ierr, "cannot receive tag values for bilinear map" )

        ierr = iMOAB_FreeSenderBuffers( cplAtmPID, atmocnid )
        call errorout( ierr, "cannot free buffers" )
    endif
   
   
   if (cplComm .ne. MPI_COMM_NULL) then

      ! We have the remapping weights now. Let us apply the weights onto the tag we defined
      ! on the source mesh and get the projection on the target mesh
      ierr = iMOAB_ApplyScalarProjectionWeights(cplAtmOcnPID, filter_type, weights_identifier1, &
                                                fields, &
                                                projectedFields)
      call errorout(ierr, 'failed to compute first order projection weight application')

      ! We have the remapping weights now. Let us apply the weights onto the tag we defined
      ! on the source mesh and get the projection on the target mesh
      ierr = iMOAB_ApplyScalarProjectionWeights(cplAtmOcnPID, filter_type, "bilinear"//C_NULL_CHAR, &
                                                fields, &
                                                projectedFieldsBilin)
      call errorout(ierr, 'failed to compute bilinear projection weight application')


      ! We have the remapping weights now. Let us apply the weights onto the tag we defined
      ! on the source mesh and get the projection on the target mesh
      ierr = iMOAB_ApplyScalarProjectionWeights(cplAtmOcnPID, filter_type, "secondorder"//C_NULL_CHAR, &
                                                fields, &
                                                projectedFieldsSecond)
      call errorout(ierr, 'failed to compute second order projection weight application')

      ! We have the remapping weights now. Let us apply the weights onto the tag we defined
      ! on the source mesh and get the projection on the target mesh
      filter_type = 1 ! global CAAS operator application
      ierr = iMOAB_ApplyScalarProjectionWeights(cplAtmOcnPID, filter_type, "secondorder"//C_NULL_CHAR, &
                                                fields, &
                                                projectedFieldsCAAS)
      call errorout(ierr, 'failed to compute second order with CAAS projection weight application')

      write(nproc,"(I0.2)")num_procs !
      outputFileOcn = "OcnOnCplF_n"//trim(nproc)//".h5m"//C_NULL_CHAR
      fileWriteOptions = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
      ierr = iMOAB_WriteMesh(cplOcnPID, outputFileOcn, fileWriteOptions)
      call errorout(ierr, 'could not write OcnOnCpl.h5m to disk')

   end if
   ! send the projected tag back to ocean pes, with send/receive tag
   ! first makje sure the tags are defined, otherwise they cannot be received
   if (ocnComm .ne. MPI_COMM_NULL) then

      ierr = iMOAB_DefineTagStorage(cmpOcnPID, transferFields, tagTypes(2), ocnCompNDoFs, tagIndexIn2)
      call errorout(ierr, 'failed to define the field tag for receiving back the tag a2oTbot_proj,  on ocn pes')
      ierr = iMOAB_DefineTagStorage(cmpOcnPID, "GLOBAL_ID"//C_NULL_CHAR, 0, 1, tagIndexIn2)
      call errorout(ierr, 'failed to define the field tag for GLOBAL_ID,  on ocn pes')

   end if

   !  send the tag to ocean pes, from ocean mesh on coupler pes
   !  from couComm, using common joint comm ocn_coupler
   !  as always, use nonblocking sends
   !  original graph (context is -1_
   if (cplComm .ne. MPI_COMM_NULL) then
      context_id = cmpocn
      ierr = iMOAB_SendElementTag(cplOcnPID, transferFields, ocnCouComm, context_id)
      call errorout(ierr, 'cannot send tag values back to ocean pes')
   end if

   if (ocnComm .ne. MPI_COMM_NULL) then
      context_id = cplocn
      ierr = iMOAB_ReceiveElementTag(cmpOcnPID, transferFields, ocnCouComm, context_id)
      call errorout(ierr, 'cannot receive tag values from ocean mesh on coupler pes')
   end if

   if (cplComm .ne. MPI_COMM_NULL) then
      context_id = cmpocn
      ierr = iMOAB_FreeSenderBuffers(cplOcnPID, context_id)
      call errorout(ierr, 'cannot free sender buffers on coupler')
   end if

   if (ocnComm .ne. MPI_COMM_NULL) then

      outputFileOcn = "OcnWithProjF.h5m"//C_NULL_CHAR
      fileWriteOptions = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
      if (my_id .eq. 0) then
         print *, ' Writing ocean mesh file with projected solution to disk: ', outputFileOcn
      end if
      ierr = iMOAB_WriteMesh(cmpOcnPID, outputFileOcn, fileWriteOptions)
      call errorout(ierr, 'could not write OcnWithProjF.h5m to disk')
      ! check baseline for some variables
      !  Each process in the communicator will have access to a local mesh instance, which
      !  will contain the original cells in the local partition and ghost entities. Number of
      !  vertices, primary cells, visible blocks, number of sidesets and nodesets boundary
      !  conditions will be returned in numProcesses 3 arrays, for local, ghost and total
      !  numbers.

      ierr = iMOAB_GetMeshInfo(cmpOcnPID, nverts, nelem, nblocks, nsbc, ndbc)
      call errorout(ierr, 'failed to get num primary elems')
      storLeng = nelem(3) ! 1 tag for now
      allocate (vals(storLeng))
      allocate (gids(storLeng))
      eetype = 1 ! double type

      eps = 1.e-9
      eetype = 1 ! cell type, not vertex
      ierr         = iMOAB_GetIntTagStorage( cmpOcnPID, "GLOBAL_ID"//C_NULL_CHAR, storLeng, eetype, gids );
      call errorout(ierr, 'failed to get gids')
      ierr         = iMOAB_GetDoubleTagStorage( cmpOcnPID, "Sa_pbot_bilin_proj"//C_NULL_CHAR, storLeng, eetype, vals );
      call errorout(ierr, 'failed to get pbots bilinear')
      
      call check_baseline(base_file3, storLeng, gids, vals, eps, my_id, ierr)
      call errorout(ierr, 'failed to check bilinear values')
      if (ierr .eq. 0 .and. my_id .eq. 0 ) print *, 'checked Sa_pbot_bilin_proj values agains baseline'
      ierr         = iMOAB_GetDoubleTagStorage( cmpOcnPID, "Sa_pbot_o2_proj"//C_NULL_CHAR, storLeng, eetype, vals )
      call errorout(ierr, 'failed to get pbot order 2')
      
      call check_baseline(base_file4, storLeng, gids, vals, eps, my_id, ierr)
      call errorout(ierr, 'failed to check higher order values')
      if (ierr .eq. 0 .and. my_id .eq. 0 ) print *, 'checked Sa_pbot_o2_proj values against baseline'

      ierr         = iMOAB_GetDoubleTagStorage( cmpOcnPID, "Sa_pbot_o2_caas_proj"//C_NULL_CHAR, storLeng, eetype, vals )
      call errorout(ierr, 'failed to get pbot order 2 caas')
      
      call check_baseline(base_file5, storLeng, gids, vals, eps, my_id, ierr)
      call errorout(ierr, 'failed to check higher order values caas')
      if (ierr .eq. 0 .and. my_id .eq. 0 ) print *, 'checked Sa_pbot_o2_caas_proj values against baseline'
      
   end if

   ! free up resources
   ierr = iMOAB_DeregisterApplication(cplAtmOcnPID)
   call errorout(ierr, 'could not de-register OCN component')

   ierr = iMOAB_DeregisterApplication(cplOcnPID)
   call errorout(ierr, 'could not de-register OCN component')
   ierr = iMOAB_DeregisterApplication(cmpOcnPID)
   call errorout(ierr, 'could not de-register OCN component')
   ierr = iMOAB_DeregisterApplication(cplAtmPID)
   call errorout(ierr, 'could not de-register OCN component')
   ierr = iMOAB_DeregisterApplication(cmpAtmPID)
   call errorout(ierr, 'could not de-register OCN component')

   ! Free all MPI datastructures
   call MPI_Comm_free(atmComm, ierr)
   call MPI_Group_free(atmGroup, ierr)
   call MPI_Comm_free(ocnComm, ierr)
   call MPI_Group_free(ocnGroup, ierr)
   call MPI_Comm_free(cplComm, ierr)
   call MPI_Group_free(cplGroup, ierr)

   call MPI_Comm_free(atmCouComm, ierr)
   call MPI_Comm_free(ocnCouComm, ierr)
   call MPI_Comm_free(global_comm, ierr)
   call MPI_Group_free(jgroup, ierr)
   call MPI_Finalize(ierr)

end program imoab_coupler_fortran
