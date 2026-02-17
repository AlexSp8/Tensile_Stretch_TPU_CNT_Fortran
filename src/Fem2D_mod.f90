
module OMP_PARALLEL

    implicit none

    integer :: NTHREADS = 20

end module OMP_PARALLEL

! ----------------------------------------------------------------------
! ======================================================================
! ----------------------------------------------------------------------

module PHYSICAL_MODULE

    implicit none

    real(8), parameter :: PI = 4.0d0*atan(1.0d0), KbN = 1.380649d-23 !in m^2 kg s^-2 K^-1
    real(8), parameter :: TempN = 25.0d0+273.15d0, PhiN = 0.04d0 !filling ratio
    real(8), parameter :: LpN = 4.0d-3, DpN = 8.0d-6 !mm
    real(8), parameter :: repN = LpN/DpN, ThetaN = (repN**2-1.0d0)/(repN**2+1.0d0)
    !UCD: +1.0d0, LCD: -1.0d0 ; Corotational: 0.0d0
    real(8), parameter :: a_UCD = 1.0d0
    !'LPTT', 'ePTT', 'Giesekus', 'FENE-CR', 'FENE-P'
    character(*), parameter :: VE_model = 'Giesekus'
    !Neo-Hookean, Mooney, Gent-Thomas, Gent-I, Gent-II
    character(*), parameter :: Elastic_model = 'Gent-I'
    !Neat, CNT, CB, ABS
    character(*), parameter :: Specimen_material = 'Neat'
    !Linear, Exponential
    character(*), parameter :: Stretch_type = 'Linear'

    integer, parameter :: Ndim = 3, Nparam = 6!5

    logical, parameter :: Read_Solution = .false.

    real(8) :: PaN         !PRESSURE AIR VISCOSITY
    real(8) :: SvN         !SOLVENT VISCOSITY
    real(8) :: UlN, EdN    !LID VELOCITY, E DOT
    real(8) :: X0, Y0, Z0, R0      !FILAMENT SIZE

    real(8), dimension(Nparam), target :: Xparam, Xparam0
    character(20), dimension(Nparam) :: Param_Names = ['TdN', 'GveN', 'GeN', 'JmN', 'EpN', 'TaN']
    integer, parameter :: Nfactors = 1
    real(8), dimension(Nfactors), parameter :: Parametric_Factors = 1.0d0
    ! integer, parameter :: Nfactors = 5
    ! real(8), dimension(Nfactors), parameter :: &
    !     Parametric_Factors = [1.0d0, 1.0d1, 1.0d2, 1.0d-1, 1.0d-2]

    type ParameterPointer
        character(20) :: name
        real(8), pointer :: p
    end type ParameterPointer

    type(ParameterPointer) :: Param

    contains

    subroutine setPhysicalParameters()

        ! use TIME_INTEGRATION_MODULE, only: dt
        use Tools_mod, only: getParameterID

        implicit none

        real(8) :: HpN, DeN, BvN
        real(8) :: GveN        !VISCOELASTIC MODULUS
        real(8) :: GeN, JmN    !ELASTIC MODULUS
        real(8) :: TdN         !REPTATION TIME
        real(8) :: EpN         !PTT PARAMETER
        real(8) :: TaN         !PARTICLE ROTATION TIME
        integer :: j

        X0 = 2.0d0/2.0d0
        Y0 = 25.0d0/2.0d0
        Z0 = 115.0d0

        PaN = 0.0d0
        SvN = 0.0d0
        ! dt = 0.01d0

        !New
        UlN = 10.0d0/60.0d0
        EdN = UlN/Z0

        !(s, mm, MPa)
        select case (Specimen_material)
        case ('Neat')
            !TdN, GveN, GeN, JmN, EpN, ThetaN, TaN
            Xparam = [2.41509412d2, 2.81306868d1, 3.66242937d0, 9.10132640d0, 0.0d0, 1.0d0] !black
            ! Xparam = [2.415185d2, 2.539215d1, 4.077423d0, 9.948906d0, 0.0d0, 1.0d0] !blue
        case ('CNT')
            !TdN, GveN, GeN, JmN, EpN, ThetaN, TaN
            ! Xparam = [2.333901d2, 2.213947d1, 6.561987d0, 4.984662d1, 0.0d0, 1.0d0] !red
            Xparam = [2.33385098d2, 2.23385289d1, 6.76035051d0, 5.05890544d1, 0.0d0, 1.0d0] !blue
            ! Xparam = [2.333858d2, 2.115874d1, 6.959257d0, 5.057385d1, 0.0d0, 1.0d0] !black
        case ('CB')
            !TdN, GveN, GeN, JmN, EpN, ThetaN, TaN
            Xparam = [2.331926d2, 2.717913d1, 4.0d0, 5.202224d1, 0.0d0, 1.0d0]
            ! Xparam = [2.326585d2, 2.685822d1, 3.66242937d0, 5.578774d1, 0.0d0, 1.0d0]
            ! Xparam = [2.325985d2, 2.773392d1, 3.66242937d0, 5.605272d1, 0.0d0, 1.0d0]
        case ('ABS')
            EdN = 5.0d0/60.0d0
            UlN = EdN*Z0
            ! Xparam = [1.008056d2, 2.641166d2, 1.167850d3, 1.0d2, 0.0d0, 1.0d0] !cold
            Xparam = [1.5d2, 2.641166d2, 1.167850d3, 1.0d2, 0.0d0, 1.0d0] !cold new
            ! Xparam = [1.000759d2, 2.510482d2, 1.151651d3, 1.0d2, 0.0d0, 1.0d0] !room
            ! Xparam = [1.005d2, 2.6d2, 1.15d3, 1.0d2, 0.0d0, 1.0d0] !room new
            ! Xparam = [1.0d2, 9.86d1, 9.95d2, 1.0d2, 0.0d0, 1.0d0] !hot

            ! Xparam = [1.0d0, 1.0d0, 2.085d3, 1.0d2, 0.0d0, 1.0d0] !cold
            ! Xparam = [1.0d0, 1.0d0, 2.07d3, 1.0d2, 0.0d0, 1.0d0] !room
            ! Xparam = [1.0d0, 1.0d0, 1.395d3, 1.0d2, 0.0d0, 1.0d0] !hot
        case default
            write(*,'(a)') 'Wrong Specimen_material in setPhysicalParameters!'
            stop
        end select
        ! Xparam = 10.0d0
        ! ! j = getParameterID('JmN')
        ! ! Xparam(j) = 5.0d0
        ! j = getParameterID('EpN')
        ! Xparam(j) = 0.0d0

        Xparam0(:) = Xparam(:)

        !--------------------------------------------------------
        ! !Varchanis: (s, mm, mPa)
        ! R0 = 1.0d0
        ! Z0 = 3.0d0*R0

        ! EdN = 1.0d0
        ! UlN = EdN*R0

        ! SvN = 1.0d0
        ! BvN = 0.262d0
        ! HpN = (1.0d0-BvN)*(SvN/BvN)
        ! DeN = 2.0d0
        ! !TdN, GveN, GeN, JmN, EpN
        ! Xparam = [DeN/EdN, HpN/TdN, 0.2d0, 3.0d0, 0.32d0]

        !--------------------------------------------------------
        ! !Stephanou: (s, mm, MPa) uniaxial elongation
        ! UlN = 2.0d0
        ! ! UlN = 10.0d0/60.0d0
        ! EdN = UlN/12.0d0
        ! DeN = 0.45d0
        ! !TdN, GveN, GeN, JmN, EpN
        ! Xparam = [DeN/EdN, 11.0d0, 0.2d0, 3.0d0, 0.0d0]

        TdN = Xparam(getParameterID('TdN'))
        GveN = Xparam(getParameterID('GveN'))
        GeN = Xparam(getParameterID('GeN'))
        JmN = Xparam(getParameterID('JmN'))
        EpN = Xparam(getParameterID('EpN'))
        TaN = Xparam(getParameterID('TaN'))

        write(*,'(2a)') 'Material: ', Specimen_material
        write(*,'(2a)') 'VE model: ', VE_model
        write(*,'(2a)') 'Elastic model: ', Elastic_model
        write(*,*)
        write(*,'(*(a,es12.4))') 'X0 = ', X0, ', Y0 = ', Y0, ', Z0 = ', Z0, ', EdN = ', EdN
        ! write(*,'(*(a,es12.4))') 'R0 = ', R0, ', Z0 = ', Z0, ', EdN = ', EdN, '
        write(*,'(*(a,es12.4))') 'SvN = ', SvN!, ', HpN = ', HpN, ', BvN = ', BvN
        write(*,'(*(a,es12.4))') 'TdN = ', TdN, ', GveN = ', GveN, ', EpN = ', EpN
        write(*,'(*(a,es12.4))') 'GeN = ', GeN, ', JmN = ', JmN
        write(*,'(*(a,es12.4))') 'repN = ', repN, ', ThetaN = ', ThetaN, ', TaN = ', TaN
        write(*,'(a)') '----------------------------------------------------------------------'

        open(15,file='Parameters.dat',position='append')
        write(15,'(*(es14.6))') Xparam, 1.0d0, 1.0d0
        close(15,status='keep')

        Param%name = ''
        Param%p => NULL()

    end subroutine setPhysicalParameters

    ! ----------------------------------------------------------------------

    real(8) function STRAIN (t)

        implicit none

        real(8), intent(in) :: t

        select case (Stretch_type)

        case ('Linear')
            STRAIN = Z0 + UlN*t
        case ('Exponential')
            STRAIN = Z0*exp(EdN*t)
        case default
            write(*,'(a)') 'Wrong Stretch_type in STRAIN!'
            stop
        end select

    end function STRAIN

end module PHYSICAL_MODULE

! ----------------------------------------------------------------------
! ======================================================================
! ----------------------------------------------------------------------

module TIME_INTEGRATION_MODULE

    implicit none

    logical, parameter :: adjust_dt = .false.
    real(8), parameter :: Initial_time = 0.0d0, Final_time = 1.0d3!3.0d-1
    real(8), parameter :: dt_min = 1.0d-5, dt_max = 1.0d0/2, Eps_pr = 1.0d-3
    ! real(8), parameter :: dt_min = 1.0d-5, dt_max = 1.0d-3, Eps_pr = 1.0d-3
    integer :: Increment
    real(8) :: Total_time, dt, dto, dtb

end module TIME_INTEGRATION_MODULE

! ----------------------------------------------------------------------
! ======================================================================
! ----------------------------------------------------------------------

module ELEMENTS_MODULE

    implicit none

    integer, parameter:: NBF_1d = 2      ! NUMBER OF BASIS FUNCTIONS PER 1d ELEMENT
    integer, parameter:: NBF_2d = 3      ! NUMBER OF BASIS FUNCTIONS PER 2d ELEMENT
    integer, parameter:: NBF_3d = 4      ! NUMBER OF BASIS FUNCTIONS PER 3d ELEMENT
    integer, parameter:: Ned_3d = 6      ! NUMBER OF EDGES PER 3d ELEMENT
    integer, parameter:: Nfc_3d = 4      ! NUMBER OF EDGES PER 3d ELEMENT
    integer, parameter:: NEQ_f = 4       ! NUMBER OF PDEs SYSTEM TO SOLVE FOR FLOW
    integer, parameter:: NEQ_s = 3*6     ! NUMBER OF PDEs SYSTEM TO SOLVE FOR STRESSES

    integer :: NEL_1d, NEL_2d, NEL_3d, NEL, NEL_fc
    integer :: NODTOL, NODTOL_fc
    integer :: NUNKNOWNS_f

    contains

    subroutine DISCRETIZATION_DATA()

        implicit none

        integer :: inod, iel, dum, index_el, idim

        open(30,file='MESH.DTA')

        read(30,*) NODTOL, NEL

        do inod = 1, NODTOL
            read(30,*) dum
        end do

        NEL_1d = 0 ; NEL_2d = 0 ; NEL_3d = 0
        do iel = 1, NEL

            read(30,*) dum, index_el

            idim = index_el/100

            if (idim == 1) then
                NEL_1d = NEL_1d + 1
            else if (idim == 2) then
                NEL_2d = NEL_2d + 1
            else if (idim == 3) then
                NEL_3d = NEL_3d + 1
            end if

        end do

        close(30)

        NUNKNOWNS_f = NODTOL*NEQ_f ! TOTAL NUMBER OF UNKNOWNS

        open(30,FILE='FACE.DTA')

        read(30,*) NODTOL_fc, NEL_fc

        close(30)

    end subroutine DISCRETIZATION_DATA

end module ELEMENTS_MODULE

! ----------------------------------------------------------------------
! ======================================================================
! ----------------------------------------------------------------------

module FLOW_ARRAYS_MODULE

    implicit none

    real(8), allocatable, dimension(:) :: B_f, S_f
    logical, allocatable, dimension(:) :: Bc_f, Bce_f

    contains

    subroutine ALLOCATE_FLOW_ARRAYS(L,NUNKNOWNS_f_)

        implicit none

        logical, intent(in) :: L
        integer, intent(in) :: NUNKNOWNS_f_

        if (L) then
            allocate( B_f(NUNKNOWNS_f_)  )  ; B_f   = 0.0d0
            allocate( S_f(NUNKNOWNS_f_)  )  ; S_f   = 0.0d0
            allocate( Bc_f(NUNKNOWNS_f_)  ) ; Bc_f  = .false.
            allocate( Bce_f(NUNKNOWNS_f_) ) ; Bce_f = .false.
        else 
            deallocate( B_f   )
            deallocate( S_f   )
            deallocate( Bc_f  )
            deallocate( Bce_f )
        end if

    end subroutine ALLOCATE_FLOW_ARRAYS

end module FLOW_ARRAYS_MODULE

! ----------------------------------------------------------------------
! ======================================================================
! ----------------------------------------------------------------------

module GLOBAL_ARRAYS_MODULE

    implicit none

    real(8), allocatable, dimension(:,:) :: TL, TLo, TLb, TLp
    real(8), allocatable, dimension(:,:) :: TD, TDo, TDb, TDp, TDg

    contains

    subroutine ALLOCATE_GLOBAL_ARRAYS(L,NODTOL_,NEL_3d_,NEQ_f_,NEQ_s_)

        implicit none

        logical, intent(in) :: L
        integer, intent(in) :: NODTOL_, NEQ_f_, NEQ_s_, NEL_3d_

        if (L) then
            allocate( TL(NODTOL_,NEQ_f_) )  ; TL  = 0.0d0
            allocate( TLo(NODTOL_,NEQ_f_) ) ; TLo = 0.0d0
            allocate( TLb(NODTOL_,NEQ_f_) ) ; TLb = 0.0d0
            allocate( TLp(NODTOL_,NEQ_f_) ) ; TLp = 0.0d0
            allocate( TD(NEL_3d_,NEQ_s_) )  ; TD  = 0.0d0
            allocate( TDo(NEL_3d_,NEQ_s_) ) ; TDo = 0.0d0
            allocate( TDb(NEL_3d_,NEQ_s_) ) ; TDb = 0.0d0
            allocate( TDp(NEL_3d_,NEQ_s_) ) ; TDp = 0.0d0
            allocate( TDg(NODTOL_,NEQ_s_) ) ; TDg = 0.0d0
        else 
            deallocate( TL  )
            deallocate( TLo )
            deallocate( TLb )
            deallocate( TLp )
            deallocate( TD  )
            deallocate( TDo )
            deallocate( TDb )
            deallocate( TDp )
            deallocate( TDg )
        end if

    end subroutine ALLOCATE_GLOBAL_ARRAYS

end module GLOBAL_ARRAYS_MODULE

! ----------------------------------------------------------------------
! ======================================================================
! ----------------------------------------------------------------------

module NRAPSHON_MODULE

    implicit none

    integer, parameter :: NITER = 50
    real(8), parameter :: ERROR_NR = 1.0d-8, ERROR_STR = 1.0d-10, MNR_eps = 1.0d-1

    integer :: iter0
    real(8) :: Cor_norm0

end module NRAPSHON_MODULE

! ----------------------------------------------------------------------
! ======================================================================
! ----------------------------------------------------------------------

module ENUMERATION_MODULE

    implicit none

    !BULK
    integer, allocatable, dimension(:,:) :: NM_MESH       ! FOR TETRAHEDRAL ELEMENTS
    integer, allocatable, dimension(:,:) :: NM_f          ! FOR TETRAHEDRAL ELEMENTS

    integer, allocatable, dimension(:,:,:) :: NME_MESH    ! FOR TETRAHEDRAL ELEMENTS
    integer, allocatable, dimension(:,:,:) :: NMF_MESH    ! FOR TETRAHEDRAL ELEMENTS

    !SURFACE
    integer, allocatable, dimension(:,:) :: NM_SURF     ! FOR TRIANGULAR ELEMENTS

    integer :: IDEsuN, IDNsuN, maxIDNsuN
    integer, allocatable, dimension(:) :: EsuN1, NsuN1
    integer, allocatable, dimension(:) :: EsuN2, NsuN2

    contains

    subroutine ALLOCATE_ENUMERATION_INDICES(L,NEL_2d_,NEL_3d_,NBF_1d_,NBF_2d_,NBF_3d_, NED_3d_, NFC_3d_)

        implicit none

        logical, intent(in) :: L
        integer, intent(in) :: NEL_2d_,NEL_3d_,NBF_1d_,NBF_2d_,NBF_3d_,NED_3d_, NFC_3d_

        if (L) then
            allocate(   NM_MESH(NEL_3d_, NBF_3d_)          ) ; NM_MESH  = 0
            allocate(      NM_f(NEL_3d_, NBF_3d_)          ) ; NM_f     = 0
            allocate(  NME_MESH(NEL_3d_, NED_3d_, NBF_1d_) ) ; NME_MESH = 0
            allocate(  NMF_MESH(NEL_3d_, NFC_3d_, NBF_2d_) ) ; NMF_MESH = 0
            allocate(   NM_SURF(NEL_2d_, NBF_2d_)          ) ; NM_SURF  = 0
        else 
            deallocate( NM_MESH   )
            deallocate( NM_f      )
            deallocate( NME_MESH  )
            deallocate( NMF_MESH  )
            deallocate( NM_SURF   )
        end if

    end subroutine ALLOCATE_ENUMERATION_INDICES

    !----------------------------------------------------------------------- 

    subroutine NNUM_ELEMENT()

        use ELEMENTS_MODULE, only: NEL_1d, NEL_2d, NEL_3d, NBF_3d, NODTOL, NBF_2d

        implicit none

        integer :: inod, iel, dum1, dum2

        open(30,FILE='MESH.DTA')

        read(30,*) dum1, dum2

        do inod = 1, NODTOL
            read(30,*) dum1
        end do

        do iel = 1, NEL_1d
            read(30,*) dum1
        end do

        do iel = 1, NEL_2d
            read(30,*) dum1, dum2, (NM_SURF(iel,inod),inod=1,NBF_2d)
        end do

        do iel = 1, NEL_3d

            read(30,*) dum1, dum2, (NM_MESH(iel,inod),inod=1,NBF_3d)

            NME_MESH(iel,1,:) = [NM_MESH(iel,1), NM_MESH(iel,2)]
            NME_MESH(iel,2,:) = [NM_MESH(iel,1), NM_MESH(iel,3)]
            NME_MESH(iel,3,:) = [NM_MESH(iel,1), NM_MESH(iel,4)]
            NME_MESH(iel,4,:) = [NM_MESH(iel,2), NM_MESH(iel,3)]
            NME_MESH(iel,5,:) = [NM_MESH(iel,3), NM_MESH(iel,4)]
            NME_MESH(iel,6,:) = [NM_MESH(iel,2), NM_MESH(iel,4)]

            NMF_MESH(iel,1,:) = [NM_MESH(iel,1), NM_MESH(iel,3), NM_MESH(iel,4)]
            NMF_MESH(iel,2,:) = [NM_MESH(iel,1), NM_MESH(iel,4), NM_MESH(iel,2)]
            NMF_MESH(iel,3,:) = [NM_MESH(iel,1), NM_MESH(iel,2), NM_MESH(iel,3)]
            NMF_MESH(iel,4,:) = [NM_MESH(iel,2), NM_MESH(iel,3), NM_MESH(iel,4)]

        end do

        close(30)

    end subroutine NNUM_ELEMENT

    !-----------------------------------------------------------------------

    integer function GNTR(GNOD)

        use ELEMENTS_MODULE, only: NEQ_f

        implicit none

        integer, intent(in) :: GNOD

        GNTR = (GNOD-1)*NEQ_f + 1

    end function GNTR

    !-----------------------------------------------------------------------

    subroutine NNUM_f()

        use ELEMENTS_MODULE, only: NEL_3d, NBF_3d

        implicit none

        integer :: iel, inod, II

        do iel = 1, NEL_3d
            do inod = 1, NBF_3d
                II = NM_MESH(iel,inod)
                NM_f(iel,inod) = GNTR(II)
            end do
        end do

    end subroutine NNUM_f

    !----------------------------------------------------------------------- 

    subroutine SURROUNDING_NUMBERING

        use ELEMENTS_MODULE, only: NBF_3d, NODTOL, NEL_3d

        implicit none

        integer :: iel, inod, JNOD, KNOD
        integer :: I, J, K, L, II, JJ, KK, LL, IM
        integer, dimension(NBF_3d) :: ISID, JSID 
        integer, allocatable, dimension(:) :: LPN
        integer, allocatable, dimension(:) :: ITEMP
        !-------------------------------------------------------------
        !    ELEMENTS SURROUNDING NODES
        !-------------------------------------------------------------
        !    EsuN1 --> ARRAY HOLDING ELEMENTS THAT SURROUND EACH NODE
        !    EsuN2 --> INDEXING ARRAY
        !    ALLOCATE ARRAY EsuN2 & LPN
        allocate( EsuN2(NODTOL+1) )
        allocate( LPN(NODTOL) )
        allocate( ITEMP(NODTOL) )

        !    COUNT NUMBER OF ELEMENTS CONNECTED TO EACH NODE
        EsuN2 = 0

        do iel  = 1, NEL_3d
            do inod = 1, NBF_3d
                II = NM_MESH(iel,inod)
                EsuN2(II) = EsuN2(II) + 1
            end do
        end do

        !    RESHUFFLING
        do II = 2, NODTOL+1
            EsuN2(II) = EsuN2(II) + EsuN2(II-1)
        end do

        !    DEFINE IDEsuN
        IDEsuN = EsuN2(NODTOL+1)

        !    RESHUFFLING
        do II = NODTOL+1, 2, -1
            EsuN2(II) = EsuN2(II-1)
        end do
        EsuN2(1) = 0

        !    ALLOCATE ARRAY EsuN1
        allocate( EsuN1(IDEsuN) )

        EsuN1 = 0

        !    STORE THE ELEMENTS IN ARRAY EsuN1
        do iel  = 1, NEL_3d
            do inod = 1, NBF_3d
                II = NM_MESH(iel,inod)
                JJ = EsuN2(II)
                LOOP_IN: do 
                    JJ = JJ + 1
                    if(EsuN1(JJ)==0)then
                        EsuN1(JJ) = iel
                        EXIT LOOP_IN
                    end if 
                end do LOOP_IN
            end do
        end do

        !-------------------------------------------------------------
        !    NODES SURROUNDING NODES
        !-------------------------------------------------------------
        !    NsuN1 --> ARRAY HOLDING NODES THAT SURROUND EACH NODE
        !    NsuN2 --> INDEXING ARRAY

        !    DEFINE maxIDNsuN
        maxIDNsuN = (NBF_3d-1)*IDEsuN

        !    ALLOCATE ARRAY NsuN2
        allocate( NsuN2(NODTOL+1) )
        allocate( NsuN1(maxIDNsuN) )

        !    INITIALIZATION
        NsuN1 = 0
        NsuN2 = 0

        !    ASSEMBLE EsuN, NsuN
        JJ = 0
        LOOP_1: do II = 1, NODTOL
            LPN = 0
            LOOP_2: do KK = EsuN2(II)+1, EsuN2(II+1)
                iel  = EsuN1(KK)
                LOOP_3: do inod = 1, NBF_3d
                    JNOD = NM_MESH(iel,inod)
                    KNOD = JNOD
                    if(LPN(JNOD)==0)then
                        JJ = JJ + 1
                        NsuN1(JJ) = KNOD 
                        LPN(JNOD)   = inod
                    end if
                end do LOOP_3

            end do LOOP_2
            NsuN2(II+1) = JJ
        end do LOOP_1

        IDNsuN = NsuN2(NODTOL+1)

        !    RESHUFFLING
        LOOPNODES: do I = 1, NODTOL

            L = 0
            do J = NsuN2(I)+1, NsuN2(I+1)
                L = L + 1
                ITEMP(L) = NsuN1(J)
            end do

            LL = 0
            do J = NsuN2(I)+1, NsuN2(I+1)
                LL = LL + 1

                JNOD = 0
                JNOD = minval(ITEMP(1:L), MASK= ITEMP(1:L) > 0)

                do KK = 1, L
                    if(JNOD == ITEMP(KK))then
                        ITEMP(KK) = 0
                    end if
                end do

                NsuN1(J) = JNOD
            end do

        end do LOOPNODES

    end subroutine SURROUNDING_NUMBERING

end module ENUMERATION_MODULE

! ----------------------------------------------------------------------
! ======================================================================
! ----------------------------------------------------------------------

module CSR_STORAGE

    implicit none

    integer :: NZ_f, N_f, NRHS_f, MAXFCT_f, MNUM_f, MTYPE_f, PHASE_f, MSGLVL_f, IDUM_f
    integer :: ERROR_f = 0
    real(8) :: DDUM_f

    real(8), allocatable, dimension(:)   :: A_f         ! LINEAR SYSTEM'S MATRIX A                                        
    integer, allocatable, dimension(:)   :: IA_f        ! INDEXING OF MATRIX A 
    integer, allocatable, dimension(:)   :: CA_f        ! COLUMNS OF MATRIX A
    integer, allocatable, dimension(:)   :: RA_f        ! ROWS OF MATRIX A
    integer, allocatable, dimension(:)   :: DA_f        ! DIAGONAL OF MATRIX A
    integer,              dimension(64)  :: IPARM_f = 0
    integer,              dimension(64)  :: PT_f    = 0

    integer, allocatable, dimension(:,:) :: CSR_f

    contains

    subroutine ALLOCATE_CSR_ARRAYS(L,K,NUNKNOWNS_f_,NZ_f_)

        implicit none

        logical, intent(in) :: L
        integer, intent(in) :: K, NUNKNOWNS_f_, NZ_f_

        if (L) then
            if (K==1) then
                allocate( IA_f(NUNKNOWNS_f_+1) )  ;  IA_f = 0
            else
                allocate(  A_f(NZ_f_) )           ;  A_f  = 0
                allocate( CA_f(NZ_f_) )           ;  CA_f = 0
                allocate( RA_f(NZ_f_) )           ;  RA_f = 0
                allocate( DA_f(NUNKNOWNS_f_) )    ;  DA_f = 0
            end if
        else 
            deallocate( IA_f )
            deallocate( CA_f )
            deallocate( A_f  )
            deallocate( RA_f )
            deallocate( DA_f )
        end if

    end subroutine ALLOCATE_CSR_ARRAYS

    !-----------------------------------------------------------------------

    subroutine ALLOCATE_CSR_CONNECTIVITY(L,NBF_3d_,NEL_3d_)

        implicit none

        logical, intent(in) :: L
        integer, intent(in) :: NBF_3d_, NEL_3d_

        if (L) then
            allocate( CSR_f(NEL_3d_,NBF_3d_*NBF_3d_) )  ;  CSR_f = 0
        else 
            deallocate( CSR_f )
        end if

    end subroutine ALLOCATE_CSR_CONNECTIVITY

    !-----------------------------------------------------------------------

    subroutine CSR_MAIN(ICHOICE)

        use ELEMENTS_MODULE, only: NEQ_f, NBF_3d, NODTOL, NEL_3d, NUNKNOWNS_f
        use ENUMERATION_MODULE, only: NM_MESH, NM_f

        implicit none

        CHARACTER(LEN=2), intent(in) :: ICHOICE

        integer :: iel, I, J

        !    NUMBERING: NU ; CONNECTIVITY: CO
        SELECT CASE(ICHOICE)

        !      DEFINE NON-ZEROS, IA & JA
        CASE('NU')
            CALL CSR_NUMBERING

        !      CSR CONNECTIVITY FOR ELEMENTAL MATRIX STORAGE
        CASE('CO')
            CALL ALLOCATE_CSR_CONNECTIVITY(.true.,NBF_3d,NEL_3d)

            !      DEFINE CSR FLOW CONNECTIVITY FOR ALL ELEMENTS
            LOOP_ALL_ELEMENTS_CSR_f: do iel = 1, NEL_3d

                CALL CSR_CONNECTIVITY&
                    (iel, NEL_3d, NM_f, NBF_3d, CSR_f, NBF_3d*NBF_3d, NZ_f, IA_f, &
                    CA_f, RA_f, NUNKNOWNS_f+1)

                CALL CSR_DIAGONAL&
                    (iel, NEL_3d, NM_f, NBF_3d, NEQ_f, DA_f, NUNKNOWNS_f, RA_f,   &
                    CA_f, NZ_f, IA_f, NUNKNOWNS_f+1)

            end do LOOP_ALL_ELEMENTS_CSR_f

        end SELECT

    end subroutine CSR_MAIN

    !-----------------------------------------------------------------------

    subroutine CSR_NUMBERING

        use ELEMENTS_MODULE, only: NODTOL, NUNKNOWNS_f, NBF_3d, NEQ_f
        use ENUMERATION_MODULE, only: IDNsuN, NsuN1, NsuN2, IDNsuN, NsuN1, NsuN2, GNTR

        implicit none

        integer :: I, II, J, JJ, K, KK, L, MZ

        NZ_f = 0

        CALL ALLOCATE_CSR_ARRAYS(.true., 1, NUNKNOWNS_f, NZ_f)

        IA_f(1) = 1  ;  JJ = 0
        do I = 1, NODTOL
            do II = 1, NEQ_f
                JJ = JJ + 1
                MZ = 0
                do J = NsuN2(I)+1, NsuN2(I+1)
                    NZ_f = NZ_f + NEQ_f
                    MZ   = MZ   + NEQ_f
                end do
                IA_f(JJ+1) = IA_f(JJ) + MZ
            end do
        end do

        CALL ALLOCATE_CSR_ARRAYS(.true., 2, NUNKNOWNS_f, NZ_f)

        do I = 1, NUNKNOWNS_f
            do J = IA_f(I), IA_f(I+1)-1
                RA_f(J) = I
            end do
        end do

        KK = 0
        do I = 1, NODTOL
            do II = 1, NEQ_f
                do J = NsuN2(I)+1, NsuN2(I+1)
                    L = NsuN1(J)
                    do JJ = 1, NEQ_f
                        KK = KK + 1
                        CA_f(KK) = GNTR(L) + JJ - 1
                    end do
                end do
            end do
        end do

    end subroutine CSR_NUMBERING

    !-----------------------------------------------------------------------

    subroutine CSR_CONNECTIVITY&
        (iel_, NEL_3d_, NM_, NBF_3d_, CSR_, EJS_, NZ_, IA_, JA_, RA_, NUNP1_)

        implicit none

        integer, intent(in) :: iel_, NEL_3d_, NBF_3d_, EJS_, NZ_, NUNP1_
        integer, dimension(NEL_3d_,NBF_3d_), intent(in)  :: NM_
        integer, dimension(NZ_),            intent(in)  :: RA_
        integer, dimension(NZ_),            intent(in)  :: JA_
        integer, dimension(NUNP1_),         intent(in)  :: IA_

        integer, dimension(NEL_3d_,EJS_),    intent(out) :: CSR_

        integer :: I, J, K, L, II, JJ

        L = 0

        do I = 1, NBF_3d_
            II = NM_(iel_,I)
            do J = 1, NBF_3d_
                JJ = NM_(iel_,J)
                L = L + 1
                do K = IA_(II), IA_(II+1)-1
                    if (RA_(K)==II .and. JA_(K)==JJ) then
                        CSR_(iel_,L) = K
                    end if
                end do

            end do
        end do

    end subroutine CSR_CONNECTIVITY

    !-----------------------------------------------------------------------

    subroutine CSR_DIAGONAL&
        (iel_, NEL_3d_, NM_f_, NBF_3d_, NEQ_f, DA_, NUNKNOWNS_, RA_, JA_, NZ_, IA_, NUNP1_)

        implicit none

        integer, intent(in) :: iel_, NEL_3d_, NBF_3d_, NUNKNOWNS_, NZ_, NUNP1_
        integer, intent(in) :: NEQ_f

        integer, dimension(NZ_),            intent(in)  :: RA_
        integer, dimension(NZ_),            intent(in)  :: JA_
        integer, dimension(NUNP1_),         intent(in)  :: IA_
        integer, dimension(NUNKNOWNS_),     intent(out) :: DA_

        integer, dimension(NEL_3d_,NBF_3d_), intent(in)  :: NM_f_

        integer :: I, J, K, II, JJ, IEQ

        do I = 1, NBF_3d_
            JJ = NM_f_(iel_,I)
            do IEQ = 0, NEQ_f-1
                II = JJ + IEQ
                do K = IA_(II), IA_(II+1)-1
                    if(RA_(K)==II .and. JA_(K)==II)then
                        DA_(II) = K
                    end if
                end do
            end do
        end do

    end subroutine CSR_DIAGONAL

    !-----------------------------------------------------------------------

    subroutine INITIALIZE_PARDISO()

        use ELEMENTS_MODULE, only: NUNKNOWNS_f
        use OMP_PARALLEL,    only: NTHREADS

        implicit none

        !    INITIALIZE PARDISO VARIABLES

        !    FOR FLOW PROBLEM
        N_f =  NUNKNOWNS_f
        NRHS_f     = 1
        MAXFCT_f   = 1
        MNUM_f     = 1
        IPARM_f(3) = NTHREADS   ! SET OMP_NUM_THREADS PARAMETER
        IPARM_f(1) = 0

        !   SETUP PARDISO CONTROL PARAMETERS AND INITIALIZE THE SOLVERS

        !   FOR FLOW PROBLEM
        MTYPE_f = 1              ! REAL & STRUCTURALLY SYMMETRIC MATRIX

        !   REORDERING AND SYMBOLIC FACTORIZATION, THIS STEP ALSO ALLOCATES
        !   ALL MEMORY THAT IS NECESSARY FOR THE FACTORIZATION

        !   FOR FLOW PROBLEM

        PHASE_f  = 11            ! ONLY REORDERING AND SYMBOLIC FACTORIZATION
        MSGLVL_f = 0             ! WITHOUT STATISTICAL INFORMATION

        CALL PARDISO&
            (PT_f, MAXFCT_f, MNUM_f, MTYPE_f, PHASE_f, N_f, A_f, IA_f, CA_f,&
            IDUM_f, NRHS_f, IPARM_f, MSGLVL_f, DDUM_f, DDUM_f, ERROR_f)

        write(*,'(a)') 'REORDERING OF A_f COMPLETED ... '

        if(ERROR_f .NE. 0)then
            print*, 'THE FOLLOWING ERROR_f WAS DETECTED: ', ERROR_f
            print*, 'DURING SYMBOLIC FACTORIZATION. SUB: PARDISO.'
            stop
        end IF

        write(*,'(a,i0)') 'NUMBER OF NONZEROS IN FACTORS   = ',IPARM_f(18)
        write(*,'(a,i0)') 'NUMBER OF FACTORIZATION MFLOPS  = ',IPARM_f(19)
        write(*,'(a)') '----------------------------------------------------------------------'

    end subroutine INITIALIZE_PARDISO

end module CSR_STORAGE

! ----------------------------------------------------------------------
! ======================================================================
! ----------------------------------------------------------------------

module MESH_MODULE

    implicit none

    real(8), allocatable, dimension(:) :: Xm, Ym, Zm
    real(8), allocatable, dimension(:) :: Nodes_bnd, X0_vals, Y0_vals, Z0_vals

    integer, parameter :: NBD = 6  !NUMBER OF DEFINED BOUNDARIES

    real(8), parameter :: EPS_MESH = 1.0d-6  !MESH TOLERANCE

    contains

    subroutine ALLOCATE_MESH_ARRAYS(L,NODTOL_, NODTOL_fc_)

        implicit none

        logical, intent(in) :: L
        integer, intent(in) :: NODTOL_, NODTOL_fc_

        if (L) then
            allocate( Xm(NODTOL_) ) ; Xm = 0.0d0
            allocate( Ym(NODTOL_) ) ; Ym = 0.0d0
            allocate( Zm(NODTOL_) ) ; Zm = 0.0d0
            allocate( X0_vals(NODTOL_fc_) ) ; X0_vals = 0.0d0
            allocate( Y0_vals(NODTOL_fc_) ) ; Y0_vals = 0.0d0
            allocate( Z0_vals(NODTOL_fc_) ) ; Z0_vals = 0.0d0
            allocate( Nodes_bnd(NODTOL_fc_) ) ; Nodes_bnd = 0
        else 
            deallocate( Xm )
            deallocate( Ym )
            deallocate( Zm )
            deallocate( X0_vals )
            deallocate( Y0_vals )
            deallocate( Z0_vals )
            deallocate( Nodes_bnd )
        end if

    end subroutine ALLOCATE_MESH_ARRAYS

    !-----------------------------------------------------------------------

    subroutine MESH()

        use ELEMENTS_MODULE, only: NODTOL, NODTOL_fc

        implicit none

        integer :: inod, dum1, dum2

        open(30,file='MESH.DTA')

        read(30,*) dum1, dum2

        do inod = 1, NODTOL
            read(30,*) dum1, Xm(inod), Ym(inod), Zm(inod)
        end do

        close(30)

        open(30,file='FACE.DTA')

        read(30,*) dum1, dum2

        do inod = 1, NODTOL_fc
            read(30,*) Nodes_bnd(inod), X0_vals(inod), Y0_vals(inod), Z0_vals(inod)
        end do

        close(30)

    end subroutine MESH

end module MESH_MODULE

! ----------------------------------------------------------------------
! ======================================================================
! ----------------------------------------------------------------------

module BOUNDARY_ENUMERATION_MODULE

    use PHYSICAL_MODULE, only: X0, Y0, R0, Z0
    use ELEMENTS_MODULE, only: NUNKNOWNS_f, NEL_3d, NODTOL, NBF_3d, NFC_3d, NBF_2d
    use ENUMERATION_MODULE, only: NM_MESH, NMF_MESH
    use MESH_MODULE, only: Xm, Ym, Zm, NBD, EPS_MESH, Nodes_bnd, X0_vals, Y0_vals, Z0_vals

    implicit none

    !BOUNDARY INDICES
    TYPE INLO
        logical, dimension(NBD) :: KOE  !  KOE(IBND) = FALSE, if ELEMENT DOES NOT BELONG TO BOUNDARY # IBND
                                        !  KOE(IBND) = TRUE,  if ELEMENT BELONGS TO BOUNDARY # IBND

        integer, dimension(NBD) :: BFE  !  BFE(IBND) = 0,   DEFAULT, THE ELEMENT DOES NOT BELONG TO A BOUNDARY
                                        !  BFE(IBND) = ifc, 1<=ifc<=NFC_3d, FACE # ifc IS ATTACHED TO BOUNDARY # IBND
        integer, dimension(NBF_2d, NBD) :: nodes
    end TYPE INLO

    TYPE(INLO), allocatable, dimension(:) :: NBE

    integer, parameter :: Nbd_Symmetry = 2
    character(1), parameter, dimension(Nbd_symmetry) :: Symmetry = ['X', 'Y']
    integer, parameter, dimension(Nbd_symmetry) :: BND_Symmetry = [1, 3]

    contains

    subroutine ALLOCATE_BOUNDARY_ENUMERATION_INDICES(L,NEL_3d_)

        implicit none

        logical, intent(in) :: L
        integer, intent(in) :: NEL_3d_

        integer :: iel

        if (L) then
            allocate( NBE(NEL_3d_) )
            do iel = 1, NEL_3d
                NBE(iel)%KOE = .false.
                NBE(iel)%BFE = 0
            end do
        else 
            deallocate( NBE )
        end if

    end subroutine ALLOCATE_BOUNDARY_ENUMERATION_INDICES

    !----------------------------------------------------------------------- 

    subroutine DEFINE_BOUNDARY_NUMBERING()

        implicit none

        integer :: iel, ifc, ibnd, N1, N2, N3, inod, gnod
        logical :: check

        do ibnd = 1, size(NBE)
            NBE(ibnd)%nodes(:,:) = 0
        end do

        ! open(8,file='nodes.dat')

        LOOP_ALL_ELEMENTS_NBE: do iel = 1, NEL_3d
            LOOP_ALL_FACES_NBE: do ifc = 1, NFC_3d

                N1 = NMF_MESH(iel,ifc,1)
                N2 = NMF_MESH(iel,ifc,2)
                N3 = NMF_MESH(iel,ifc,3)

                !X=0
                ibnd = 1
                if ((abs(Xm(N1)-0.0d0)<EPS_MESH).and.&
                    (abs(Xm(N2)-0.0d0)<EPS_MESH).and.&
                    (abs(Xm(N3)-0.0d0)<EPS_MESH)) then

                    NBE(iel)%KOE(ibnd) = .true.
                    NBE(iel)%BFE(ibnd) = ifc
                    NBE(iel)%nodes(:,ibnd) = [N1,N2,N3]

                end if

                !X=X0
                ibnd = 2
                if ((abs( Xm(N1)-X0 )<EPS_MESH).and.&
                    (abs( Xm(N2)-X0 )<EPS_MESH).and.&
                    (abs( Xm(N3)-X0 )<EPS_MESH)) then

                    NBE(iel)%KOE(ibnd) = .true.
                    NBE(iel)%BFE(ibnd) = ifc
                    NBE(iel)%nodes(:,ibnd) = [N1,N2,N3]

                end if

                !Y=0
                ibnd = 3
                if ((abs(Ym(N1)-0.0d0)<EPS_MESH).and.&
                    (abs(Ym(N2)-0.0d0)<EPS_MESH).and.&
                    (abs(Ym(N3)-0.0d0)<EPS_MESH)) then

                    NBE(iel)%KOE(ibnd) = .true.
                    NBE(iel)%BFE(ibnd) = ifc
                    NBE(iel)%nodes(:,ibnd) = [N1,N2,N3]

                end if

                ! !Y=Y0
                ! ibnd = 4
                ! if ((abs( Ym(N1)-Y0 )<EPS_MESH).and.&
                !     (abs( Ym(N2)-Y0 )<EPS_MESH).and.&
                !     (abs( Ym(N3)-Y0 )<EPS_MESH)) then

                !     NBE(iel)%KOE(ibnd) = .true.
                !     NBE(iel)%BFE(ibnd) = ifc
                !     NBE(iel)%nodes(:,ibnd) = [N1,N2,N3]

                ! end if

                !Y=Yi
                ibnd = 4
                do inod = 1, size(Nodes_bnd)

                    gnod = Nodes_bnd(inod)

                    check = any(NMF_MESH(iel,ifc,:)==gnod)

                    ! check = (abs( Xm(N1)-X0_vals(inod) )<EPS_MESH).and. &
                    !         (abs( Xm(N2)-X0_vals(inod) )<EPS_MESH).and. &
                    !         (abs( Xm(N3)-X0_vals(inod) )<EPS_MESH).and. &
                    !         (abs( Ym(N1)-Y0_vals(inod) )<EPS_MESH).and. &
                    !         (abs( Ym(N2)-Y0_vals(inod) )<EPS_MESH).and. &
                    !         (abs( Ym(N3)-Y0_vals(inod) )<EPS_MESH).and. &
                    !         (abs( Zm(N1)-Z0_vals(inod) )<EPS_MESH).and. &
                    !         (abs( Zm(N2)-Z0_vals(inod) )<EPS_MESH).and. &
                    !         (abs( Zm(N3)-Z0_vals(inod) )<EPS_MESH)

                    if (check) then

                        NBE(iel)%KOE(ibnd) = .true.
                        NBE(iel)%BFE(ibnd) = ifc
                        NBE(iel)%nodes(:,ibnd) = [N1,N2,N3]

                        ! write(8,*) inod, gnod, iel+NEL_2d+NEL_1d+NODTOL+1

                        exit

                    end if

                end do

                !Z=0
                ibnd = 5
                if ((abs(Zm(N1)-0.0d0)<EPS_MESH).and.&
                    (abs(Zm(N2)-0.0d0)<EPS_MESH).and.&
                    (abs(Zm(N3)-0.0d0)<EPS_MESH)) then

                    NBE(iel)%KOE(ibnd) = .true.
                    NBE(iel)%BFE(ibnd) = ifc
                    NBE(iel)%nodes(:,ibnd) = [N1,N2,N3]

                end if

                !Z=Z0
                ibnd = 6
                if ((abs(Zm(N1)-Z0)<EPS_MESH).and.&
                    (abs(Zm(N2)-Z0)<EPS_MESH).and.&
                    (abs(Zm(N3)-Z0)<EPS_MESH)) then

                    NBE(iel)%KOE(ibnd) = .true. 
                    NBE(iel)%BFE(ibnd) = ifc
                    NBE(iel)%nodes(:,ibnd) = [N1,N2,N3]

                end if

            end do LOOP_ALL_FACES_NBE
        end do LOOP_ALL_ELEMENTS_NBE

    end subroutine DEFINE_BOUNDARY_NUMBERING

end module BOUNDARY_ENUMERATION_MODULE

! ----------------------------------------------------------------------
! ======================================================================
! ----------------------------------------------------------------------

module GAUSS_MODULE

    use ELEMENTS_MODULE, only: NED_3d, NFC_3d, NBF_3d  

    implicit none

    integer, parameter :: NGAUSS_1d = 1
    integer, parameter :: NGAUSS_2d = 3
    integer, parameter :: NGAUSS_3d = 4

    real(8), dimension(NGAUSS_1d) :: GAPT_1d
    real(8), dimension(NGAUSS_1d) :: WO_1d

    real(8), dimension(NGAUSS_2d,2) :: GAPT_2d
    real(8), dimension(NGAUSS_2d) :: WO_2d

    real(8), dimension(NGAUSS_3d,3) :: GAPT_3d
    real(8), dimension(NGAUSS_3d) :: WO_3d

    real(8), dimension(NBF_3d,NGAUSS_1d,NED_3d) :: BFN_E, DFDC_E, DFDE_E, DFDS_E
    real(8), dimension(NBF_3d,NGAUSS_2d,NFC_3d) :: BFN_F, DFDC_F, DFDE_F, DFDS_F
    real(8), dimension(NBF_3d,NGAUSS_3d) :: BFN_3d, DFDC_3d, DFDE_3d, DFDS_3d

    contains

    subroutine setGaussIntegration()

        CALL setGaussIntegration_1d()
        CALL setGaussIntegration_2d()
        CALL setGaussIntegration_3d()

    end subroutine setGaussIntegration

    !----------------------------------------------------------------------- 

    subroutine setGaussIntegration_1d

    implicit none

        integer :: i, ig

        real(8), dimension(NBF_3d) :: BFN, DFDC, DFDE, DFDS

        GAPT_1d(1) = 0.5d0

        WO_1d(1)   = 1.0d0

        ig = 0

        do i = 1, NGAUSS_1d

            ig = ig  + 1

            !EDGE 1 (C,E=0,S=0)
            CALL LINEAR_TETRAHEDRON( BFN, DFDC, DFDE, DFDS, GAPT_1d(i), 0.0d0, 0.0d0 )
            BFN_E(:,ig,1)  = BFN
            DFDC_E(:,ig,1) = DFDC
            DFDE_E(:,ig,1) = DFDE
            DFDS_E(:,ig,1) = DFDS

            !EDGE 2 (C=0,E,S=0)
            CALL LINEAR_TETRAHEDRON( BFN, DFDC, DFDE, DFDS, 0.0d0, GAPT_1d(i), 0.0d0 )
            BFN_E(:,ig,2)  = BFN
            DFDC_E(:,ig,2) = DFDC
            DFDE_E(:,ig,2) = DFDE
            DFDS_E(:,ig,2) = DFDS

            !EDGE 3 (C=0,E=0,S)
            CALL LINEAR_TETRAHEDRON( BFN, DFDC, DFDE, DFDS, 0.0d0, 0.0d0, GAPT_1d(i) )
            BFN_E(:,ig,3)  = BFN
            DFDC_E(:,ig,3) = DFDC
            DFDE_E(:,ig,3) = DFDE
            DFDS_E(:,ig,3) = DFDS

            !EDGE 4 (C,E=1-C,S=0)
            CALL LINEAR_TETRAHEDRON( BFN, DFDC, DFDE, DFDS, GAPT_1d(i), 1.0d0-GAPT_1d(i), 0.0d0 )
            BFN_E(:,ig,4)  = BFN
            DFDC_E(:,ig,4) = DFDC
            DFDE_E(:,ig,4) = DFDE
            DFDS_E(:,ig,4) = DFDS

            !EDGE 5 (C=0,E,S=1-E)
            CALL LINEAR_TETRAHEDRON( BFN, DFDC, DFDE, DFDS, 0.0d0, GAPT_1d(i), 1.0d0-GAPT_1d(i) )
            BFN_E(:,ig,5)  = BFN
            DFDC_E(:,ig,5) = DFDC
            DFDE_E(:,ig,5) = DFDE
            DFDS_E(:,ig,5) = DFDS

            !EDGE 6 (C,E=0,S=1-C)
            CALL LINEAR_TETRAHEDRON( BFN, DFDC, DFDE, DFDS, GAPT_1d(i), 0.0d0, 1.0d0-GAPT_1d(i) )
            BFN_E(:,ig,6)  = BFN
            DFDC_E(:,ig,6) = DFDC
            DFDE_E(:,ig,6) = DFDE
            DFDS_E(:,ig,6) = DFDS

        end do

    end subroutine setGaussIntegration_1d

    !----------------------------------------------------------------------- 

    subroutine setGaussIntegration_2d

        integer :: i, ig

        real(8), dimension(NBF_3d) :: BFN, DFDC, DFDE, DFDS

        GAPT_2d(1,1) = 1.0d0/2.0d0
        GAPT_2d(1,2) = 1.0d0/2.0d0

        GAPT_2d(2,1) = 0.0d0
        GAPT_2d(2,2) = 1.0d0/2.0d0

        GAPT_2d(3,1) = 1.0d0/2.0d0
        GAPT_2d(3,2) = 0.0d0

        WO_2d(1) = (1.0d0/3.0d0)*(1.0d0/2.0d0)
        WO_2d(2) = (1.0d0/3.0d0)*(1.0d0/2.0d0)
        WO_2d(3) = (1.0d0/3.0d0)*(1.0d0/2.0d0)

        ig = 0

        do i = 1, NGAUSS_2d

            ig = ig  + 1

            !FACE 1 (C=0)
            CALL LINEAR_TETRAHEDRON( BFN, DFDC, DFDE, DFDS, 0.0d0, GAPT_2d(i,1), GAPT_2d(i,2) )
            BFN_F(:,ig,1)  = BFN
            DFDC_F(:,ig,1) = DFDC
            DFDE_F(:,ig,1) = DFDE
            DFDS_F(:,ig,1) = DFDS

            !FACE 2 (E=0)
            CALL LINEAR_TETRAHEDRON( BFN, DFDC, DFDE, DFDS, GAPT_2d(i,1), 0.0d0, GAPT_2d(i,2) )
            BFN_F(:,ig,2)  = BFN
            DFDC_F(:,ig,2) = DFDC
            DFDE_F(:,ig,2) = DFDE
            DFDS_F(:,ig,2) = DFDS

            !FACE 3 (S=0)
            CALL LINEAR_TETRAHEDRON( BFN, DFDC, DFDE, DFDS, GAPT_2d(i,1), GAPT_2d(i,2), 0.0d0 )
            BFN_F(:,ig,3)  = BFN
            DFDC_F(:,ig,3) = DFDC
            DFDE_F(:,ig,3) = DFDE
            DFDS_F(:,ig,3) = DFDS

            !FACE 4 (S=1-C-E)
            CALL LINEAR_TETRAHEDRON( BFN, DFDC, DFDE, DFDS, GAPT_2d(i,1), GAPT_2d(i,2), 1.0d0-GAPT_2d(i,1)-GAPT_2d(i,2) )
            BFN_F(:,ig,4)  = BFN
            DFDC_F(:,ig,4) = DFDC
            DFDE_F(:,ig,4) = DFDE
            DFDS_F(:,ig,4) = DFDS

        end do

    end subroutine setGaussIntegration_2d

    !----------------------------------------------------------------------- 

    subroutine setGaussIntegration_3d

        integer :: i, ig

        real(8), dimension(NBF_3d) :: BFN, DFDC, DFDE, DFDS

        GAPT_3d(1,1) = (5.0d0-sqrt(5.0d0))/20.0d0
        GAPT_3d(1,2) = (5.0d0-sqrt(5.0d0))/20.0d0
        GAPT_3d(1,3) = (5.0d0-sqrt(5.0d0))/20.0d0

        GAPT_3d(2,1) = (5.0d0-sqrt(5.0d0))/20.0d0
        GAPT_3d(2,2) = (5.0d0-sqrt(5.0d0))/20.0d0
        GAPT_3d(2,3) = (5.0d0+3.0d0*sqrt(5.0d0))/20.0d0

        GAPT_3d(3,1) = (5.0d0-sqrt(5.0d0))/20.0d0
        GAPT_3d(3,2) = (5.0d0+3.0d0*sqrt(5.0d0))/20.0d0
        GAPT_3d(3,3) = (5.0d0-sqrt(5.0d0))/20.0d0

        GAPT_3d(4,1) = (5.0d0+3.0d0*sqrt(5.0d0))/20.0d0
        GAPT_3d(4,2) = (5.0d0-sqrt(5.0d0))/20.0d0
        GAPT_3d(4,3) = (5.0d0-sqrt(5.0d0))/20.0d0

        WO_3d(1) = (1.0d0/4.0d0)*(1.0d0/6.0d0)
        WO_3d(2) = (1.0d0/4.0d0)*(1.0d0/6.0d0)
        WO_3d(3) = (1.0d0/4.0d0)*(1.0d0/6.0d0)
        WO_3d(4) = (1.0d0/4.0d0)*(1.0d0/6.0d0)

        ig = 0

        do i = 1, NGAUSS_3d

            ig = ig  + 1

            CALL LINEAR_TETRAHEDRON( BFN, DFDC, DFDE, DFDS, GAPT_3d(i,1), GAPT_3d(i,2), GAPT_3d(i,3) )

            BFN_3d(:,ig)  = BFN
            DFDC_3d(:,ig) = DFDC ! C-> ksi
            DFDE_3d(:,ig) = DFDE ! E-> eta
            DFDS_3d(:,ig) = DFDS ! S-> zeta

        end do

    end subroutine setGaussIntegration_3d

    !----------------------------------------------------------------------- 

    subroutine LINEAR_TETRAHEDRON( BFN, DFDC, DFDE, DFDS, C, E, S )

        implicit none

        real(8), intent(in) :: C, E, S

        real(8), INTENT(inout), dimension(:) :: BFN, DFDC, DFDE, DFDS

        BFN(1) = 1.0d0-C-E-S
        BFN(2) = C
        BFN(3) = E
        BFN(4) = S

        DFDC(1) = -1.0d0
        DFDC(2) =  1.0d0
        DFDC(3) =  0.0d0
        DFDC(4) =  0.0d0

        DFDE(1) = -1.0d0
        DFDE(2) =  0.0d0
        DFDE(3) =  1.0d0
        DFDE(4) =  0.0d0

        DFDS(1) = -1.0d0
        DFDS(2) =  0.0d0
        DFDS(3) =  0.0d0
        DFDS(4) =  1.0d0

    end subroutine LINEAR_TETRAHEDRON

end module GAUSS_MODULE

! ----------------------------------------------------------------------
! ======================================================================
! ----------------------------------------------------------------------

module VISUALIZATION_MODULE

    implicit none
    
    contains

    subroutine READ_SOLUTION()

        use ELEMENTS_MODULE, only: NEL_3d, NODTOL, NBF_3d
        use ENUMERATION_MODULE, only: NM_MESH
        use GLOBAL_ARRAYS_MODULE, only: TLo, TDo, TDg

        implicit none

        integer :: i, iel, inod

        open(14, file='SOLUTION.DTA')

        read(14,*)
        read(14,*)
        read(14,*)
        read(14,*)
        read(14,*)
        read(14,*)
        read(14,*)
        read(14,*)
        read(14,*)
        read(14,*)
        read(14,*)
        read(14,*)
        read(14,*)
        read(14,*)
        read(14,*)
        read(14,*)

        read(14,*) (TLo(i,1),i=1,NODTOL)
        read(14,*) (TLo(i,2),i=1,NODTOL)
        read(14,*) (TLo(i,3),i=1,NODTOL)
        read(14,*) (TLo(i,4),i=1,NODTOL)
        read(14,*) (TDg(i,1),i=1,NODTOL)
        read(14,*) (TDg(i,2),i=1,NODTOL)
        read(14,*) (TDg(i,3),i=1,NODTOL)
        read(14,*) (TDg(i,4),i=1,NODTOL)
        read(14,*) (TDg(i,5),i=1,NODTOL)
        read(14,*) (TDg(i,6),i=1,NODTOL)

        close(14)

        !    DEFINE DISCONTINOUS STRESSES
        TDo = 0.0d0
        do iel = 1, NEL_3d
            do inod = 1, NBF_3d
                TDo(iel,:) = TDo(iel,:) + TDg(NM_MESH(iel,inod),:)
            end do
            TDo(iel,:) = TDo(iel,:)/DBLE(NBF_3d)
        end do

    end subroutine READ_SOLUTION

    !----------------------------------------------------------------------- 

    subroutine writeSTLFile( ITER )

        use ELEMENTS_MODULE, only: NBF_2d, NEL_2d
        use ENUMERATION_MODULE, only: NM_SURF
        use GLOBAL_ARRAYS_MODULE, only: TL

        implicit none

        integer, intent(in) :: ITER

        integer :: ierror, iel, gnod, inod
        CHARACTER(LEN=100) :: str, FN
        CHARACTER(LEN=:), allocatable :: STR_ITER

        write(str,'(i8)') ITER
        STR_ITER = trim(adjustl(str))

        FN = 'SURF_'//STR_ITER//'.STL'

        FN = trim(adjustl(FN))

        open(14,file=FN,  status='UNKNOWN', action='WRITE', iostat=ierror,position='REWIND')

        if ( ierror /= 0 ) then
            write(*,*)'Error:  SURF.STL'
            stop
        end if

        write(14,*)'solid Mesh_1'

        do iel = 1, NEL_2d
            write(14,*)'facet normal  0.00000000 -0.00000000 0.00000000'

            write(14,*)'  outer loop'

            do inod = 1, NBF_2d
                gnod = NM_SURF(iel,inod)

                write(14,70) TL(gnod,1), TL(gnod,2), TL(gnod,3)
            end do

            write(14,*)'  endloop'
            write(14,*)'endfacet'
        end do

        write(14,*)'endsolid Mesh_1'

        close( 14, status='keep' )

        70 format('    vertex', 2x, f32.16, 2x, f32.16, 2x, f32.16)

    end subroutine writeSTLFile

end module VISUALIZATION_MODULE
