! ----------------------------------------------------------------------
!    FEM IS THE MAIN PROGRAM FOR THE SOLUTION OF A SYSTEM OF PDEs
!    IN A CUBIC DOMAIN X,Y,Z
!    IT USES TRILINEAR BASIS FUNCTIONS & HEXAHEDRAL ELEMENTS
! ----------------------------------------------------------------------
!    PROGRAM WRITTEN BY J.A.T. & A.N.M. (1985, MIT    ) 
!    REVISED BY Y.D.                    (2013, UPATRAS)
!    REVISED BY S.E.V.                  (2018, UPATRAS)
! ----------------------------------------------------------------------

program FEM2D

    use PHYSICAL_MODULE, only: setPhysicalParameters
    use ELEMENTS_MODULE, only: NODTOL, NODTOL_fc, NEL_2d, NEL_3d, &
        NBF_1d, NBF_2d, NBF_3d, NED_3d, NFC_3d, NUNKNOWNS_f, NEQ_f, NEQ_s, &
        DISCRETIZATION_DATA
    use MESH_MODULE, only: ALLOCATE_MESH_ARRAYS, MESH
    use GAUSS_MODULE, only: setGaussIntegration
    use ENUMERATION_MODULE, only: ALLOCATE_ENUMERATION_INDICES, &
        NNUM_ELEMENT, SURROUNDING_NUMBERING, NNUM_f
    use BOUNDARY_ENUMERATION_MODULE, only: ALLOCATE_BOUNDARY_ENUMERATION_INDICES, &
        DEFINE_BOUNDARY_NUMBERING
    use CSR_STORAGE, only: NZ_f, CSR_MAIN, INITIALIZE_PARDISO, &
        ALLOCATE_CSR_ARRAYS, ALLOCATE_CSR_CONNECTIVITY
    use FLOW_ARRAYS_MODULE, only: ALLOCATE_FLOW_ARRAYS
    use GLOBAL_ARRAYS_MODULE, only: ALLOCATE_GLOBAL_ARRAYS
    use Tools_mod, only: writeElapsedTime
    use Regression_mod, only: Regression, allocateRegressionSystem
    use Continuation_mod, only: solveProblem
    use omp_lib, only: omp_get_wtime

    implicit none

    character(len=50) :: date, host
    real(8) :: tStart, tEnd

    tStart = omp_get_wtime()

    call execute_command_line('mkdir out/')
    call execute_command_line('mkdir out/Tec/')
    ! call execute_command_line('mkdir Stl')
    call execute_command_line('cp -r ConvertToBinary out/Tec/')
    call execute_command_line ('mkdir -p Graphs/')
    call execute_command_line ('mkdir -p out/Graphs_full/')
    call execute_command_line ('mkdir -p out/SOL/')

    write(*,*)
    write(*,'(a)') '----------------------------------------------------------------------'

    call fdate(date)
    write(*,'(a)') date
    write(*,*)

    call execute_command_line('hostname >> host.txt')
    open(12, file='host.txt')
    read(12,*) host
    close(12, status='delete')
    write(*,'(2a)') 'Running in: ', host
    write(*,*)

    call setPhysicalParameters()

    call setGaussIntegration()

    call DISCRETIZATION_DATA()

    call ALLOCATE_MESH_ARRAYS(.true., NODTOL, NODTOL_fc)
    call MESH()

    call ALLOCATE_ENUMERATION_INDICES(.true.,NEL_2d,NEL_3d,NBF_1d,NBF_2d,NBF_3d,NED_3d,NFC_3d)
    call NNUM_ELEMENT()
    call SURROUNDING_NUMBERING()
    call NNUM_f()
    call ALLOCATE_BOUNDARY_ENUMERATION_INDICES(.true.,NEL_3d)
    call DEFINE_BOUNDARY_NUMBERING()

    call CSR_MAIN('NU')
    call CSR_MAIN('CO')
    call INITIALIZE_PARDISO()

    call ALLOCATE_FLOW_ARRAYS( .true., NUNKNOWNS_f )
    call ALLOCATE_GLOBAL_ARRAYS( .true., NODTOL, NEL_3d, NEQ_f, NEQ_s )

    if (Regression) then
        call allocateRegressionSystem()
    end if

    tEnd = omp_get_wtime()
    call writeElapsedTime(tStart,tEnd,'Initialization time')

    call solveProblem()

    call ALLOCATE_MESH_ARRAYS(.false.,NODTOL,NODTOL_fc)
    call ALLOCATE_ENUMERATION_INDICES(.false.,NEL_2d,NEL_3d,NBF_1d,NBF_2d,NBF_3d,NED_3d,NFC_3d)
    call ALLOCATE_FLOW_ARRAYS(.false.,NUNKNOWNS_f)
    call ALLOCATE_GLOBAL_ARRAYS( .false., NODTOL, NEL_3d, NEQ_f, NEQ_s )
    call ALLOCATE_CSR_ARRAYS(.false., 1, NUNKNOWNS_f, NZ_f)
    call ALLOCATE_CSR_CONNECTIVITY(.false.,NBF_3d,NEL_3d)
    call ALLOCATE_BOUNDARY_ENUMERATION_INDICES(.false.,NEL_3d)

    tEnd = omp_get_wtime()
    call writeElapsedTime(tStart,tEnd,'Total time')

end program FEM2D
