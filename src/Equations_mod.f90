
module Equations_mod

    implicit none

    contains

    subroutine FLOW_EQUATIONS(iel, Flag_NR)

        use ELEMENTS_MODULE, only: NBF_3d, NEQ_f, NEQ_s
        use ENUMERATION_MODULE, only: NM_MESH
        use BOUNDARY_ENUMERATION_MODULE, only: NBE
        use GLOBAL_ARRAYS_MODULE, only: TL, TD

        implicit none

        integer, intent(in) :: iel
        character(*), intent(in) :: Flag_NR

        integer :: inod, gnod, ibnd

        real(8), dimension(NBF_3d,NEQ_f) :: TEMP_TL
        real(8), dimension(NBF_3d,NEQ_f) :: TEMP_RES
        real(8), dimension(NEQ_s) :: TEMP_TD

        TEMP_TL(:,:) = 0.0d0
        do inod = 1, NBF_3d
            gnod = NM_MESH(iel,inod)
            TEMP_TL(inod,:) = TL(gnod,:)
        end do

        TEMP_TD(:) = TD(iel,:)

        call bulkResidual( iel, TEMP_TL, TEMP_TD, TEMP_RES, .true. )
        if (Flag_NR == 'NRP') call bulkJacobian( iel, TEMP_TL, TEMP_TD, TEMP_RES )

        ! !Y=Y0
        ! ibnd = 4
        ! if (NBE(iel)%KOE(ibnd)) then
        !     call INTERFACE_FORCE_RESIDUAL_f( iel, NBE(iel)%BFE(ibnd), TEMP_TL, TEMP_RES, .true. )
        !     if (Flag_NR == 'NRP') call INTERFACE_FORCE_JACOBIAN_f( iel, NBE(iel)%BFE(ibnd), TEMP_TL, TEMP_RES )
        ! end if

    end subroutine FLOW_EQUATIONS

    !---------------------------------------------------------------------

    subroutine bulkResidual( iel, TEMP_TL, TEMP_TD, TEMP_RES, STORE )

        use TIME_INTEGRATION_MODULE, only: INCREMENT, dt
        use PHYSICAL_MODULE, only: SvN, UlN, Xparam, EdN
        use ELEMENTS_MODULE, only: NBF_3d, NEQ_f, NEQ_s, NUNKNOWNS_f
        use GAUSS_MODULE, only: WO_3d, NGAUSS_3d, BFN_3d, DFDC_3d, DFDE_3d, DFDS_3d
        use ENUMERATION_MODULE, only: NM_MESH, NM_f
        use GLOBAL_ARRAYS_MODULE, only: TLo, TLb, TD, TDo, TDb
        use FLOW_ARRAYS_MODULE, only: B_f
        ! use MESH_MODULE, only: Xm, Ym, Zm
        use Tools_mod, only: doubleDotProduct, determinant3D, traceTensor, getParameterID, &
            storeResidual, basis3D
        use Stress_mod, only: DOMI_STRESS, VE_tensor, Elastic_tensor, Orientation_tensor

        implicit none

        integer, intent(in) :: iel
        real(8), dimension(NBF_3d, NEQ_f), intent(in) :: TEMP_TL
        real(8), dimension(NEQ_s), intent(inout) :: TEMP_TD
        real(8), dimension(NBF_3d, NEQ_f), intent(out) :: TEMP_RES
        logical, intent(in) :: STORE

        integer, parameter, dimension(3,3) :: I_Tensor = [1,0,0,0,1,0,0,0,1]
        integer :: ig, inod, gnod, i
        real(8) :: GveN, GeN, P, Continuity, trace_B, detS
        real(8), dimension(3) :: dPdXi, Mom, div_P
        real(8), dimension(3,3) :: C_Tensor, Tve_Tensor, G_dot, GV, GVT, P_Tensor, A_Tensor, Ta_Tensor
        real(8), dimension(3,3) :: GU, GUT, F_Tensor, FT_Tensor, B_Tensor, Te_Tensor, dev_B, S_Tensor
        real(8), dimension(3,3) :: dXidci
        real(8), dimension(3) :: Xi_gs, dt_term
        real(8), dimension(NEQ_f)  :: TERM_RES
        integer, dimension(NBF_3d) :: NM
        real(8), dimension(NBF_3d) :: bfn
        real(8), dimension(NBF_3d,3) :: Xi_el, dfdci, dfdxi
        real(8) :: WET, CJAC, E_TH, Helem, Ha, tpspg, tlsic
        real(8) :: BIFN, DBIX, DBIY, DBIZ

        Xi_el(:,:) = TEMP_TL(:,1:3)
        ! Xi_el(:,1) = Xm(NM_MESH(iel,:))
        ! Xi_el(:,2) = Ym(NM_MESH(iel,:))
        ! Xi_el(:,3) = Zm(NM_MESH(iel,:))

        E_TH = 0.0d0
        do ig = 1, NGAUSS_3d

            bfn(:) = BFN_3d(:,ig)
            dfdci(:,1) = DFDC_3d(:,ig)
            dfdci(:,2) = DFDE_3d(:,ig)
            dfdci(:,3) = DFDS_3d(:,ig)

            call basis3D(Xi_el, bfn, dfdci, dfdxi, Xi_gs, dXidci, CJAC)

            WET = WO_3d(ig)*abs(CJAC)

            E_TH = E_TH + WET

        end do
        Helem = E_TH**(1.0d0/3.0d0) 
        
        call DOMI_STRESS( iel, TEMP_TL,  TEMP_TD )

        C_Tensor(1,1) = TEMP_TD(1) ; C_Tensor(1,2) = TEMP_TD(2) ; C_Tensor(1,3) = TEMP_TD(3)
        C_Tensor(2,1) = TEMP_TD(2) ; C_Tensor(2,2) = TEMP_TD(4) ; C_Tensor(2,3) = TEMP_TD(5)
        C_Tensor(3,1) = TEMP_TD(3) ; C_Tensor(3,2) = TEMP_TD(5) ; C_Tensor(3,3) = TEMP_TD(6)

        S_Tensor(1,1) = TEMP_TD(7) ; S_Tensor(1,2) = TEMP_TD(8)  ; S_Tensor(1,3) = TEMP_TD(9)
        S_Tensor(2,1) = TEMP_TD(8) ; S_Tensor(2,2) = TEMP_TD(10) ; S_Tensor(2,3) = TEMP_TD(11)
        S_Tensor(3,1) = TEMP_TD(9) ; S_Tensor(3,2) = TEMP_TD(11) ; S_Tensor(3,3) = TEMP_TD(12)

        A_Tensor(1,1) = TEMP_TD(13) ; A_Tensor(1,2) = TEMP_TD(14) ; A_Tensor(1,3) = TEMP_TD(15)
        A_Tensor(2,1) = TEMP_TD(14) ; A_Tensor(2,2) = TEMP_TD(16) ; A_Tensor(2,3) = TEMP_TD(17)
        A_Tensor(3,1) = TEMP_TD(15) ; A_Tensor(3,2) = TEMP_TD(17) ; A_Tensor(3,3) = TEMP_TD(18)

        call VE_tensor(C_Tensor, Tve_Tensor)
        call Elastic_tensor(S_Tensor, Te_Tensor)
        call Orientation_tensor(A_Tensor, Ta_Tensor)

        detS = determinant3D(S_Tensor)

        TEMP_RES(:,:) = 0.0d0
        LOOP_GAUSS: do ig = 1, NGAUSS_3d

            bfn(:) = BFN_3d(:,ig)
            dfdci(:,1) = DFDC_3d(:,ig)
            dfdci(:,2) = DFDE_3d(:,ig)
            dfdci(:,3) = DFDS_3d(:,ig)

            call basis3D(Xi_el, bfn, dfdci, dfdxi, Xi_gs, dXidci, CJAC)

            WET = WO_3d(ig)*abs(CJAC)

            GU(:,:) = 0.0d0
            P = 0.0d0 ; dPdXi(:) = 0.0d0
            GV(:,:) = 0.0d0

            do inod = 1, NBF_3d

                gnod = NM_MESH(iel,inod)

                ! GU(1,1) = GU(1,1) + TEMP_TL(inod,1)*dfdxi(inod,1)
                ! GU(2,1) = GU(2,1) + TEMP_TL(inod,1)*dfdxi(inod,2)
                ! GU(1,2) = GU(1,2) + TEMP_TL(inod,2)*dfdxi(inod,1)
                ! GU(2,2) = GU(2,2) + TEMP_TL(inod,2)*dfdxi(inod,2)

                P = P + TEMP_TL(inod,4)*bfn(inod)
                dPdXi(1) = dPdXi(1) + TEMP_TL(inod,4)*dfdxi(inod,1)
                dPdXi(2) = dPdXi(2) + TEMP_TL(inod,4)*dfdxi(inod,2)
                dPdXi(3) = dPdXi(3) + TEMP_TL(inod,4)*dfdxi(inod,3)

                if (Increment == 1) then
                    dt_term(:) = (TEMP_TL(inod,1:3)-TLo(gnod,1:3))/dt
                else
                    dt_term(:) = (3.0d0*TEMP_TL(inod,1:3)-4.0d0*TLo(gnod,1:3)+TLb(gnod,1:3))/(2.0d0*dt)
                end if
                do i = 1, 3
                    GV(:,i) = GV(:,i) + dt_term(i)*dfdxi(inod,:)
                end do

            end do

            GVT(:,:) = transpose(GV)
            G_dot(:,:) = GV(:,:) + GVT(:,:)

            GveN = Xparam(getParameterID('GveN'))
            GeN = Xparam(getParameterID('GeN'))

            ! !Solid: F Tensor inverse notation
            ! GUT(:,:) = transpose(GU)
            ! !B Tensor
            ! F_Tensor(:,:) = I_Tensor(:,:) + GUT(:,:)
            ! FT_Tensor(:,:) = transpose(F_Tensor)
            ! ! FT_Tensor(:,:) = I_Tensor(:,:) + GU(:,:)
            ! B_Tensor(:,:) = matmul( FT_Tensor, F_Tensor )
            ! trace_B = B_Tensor(1,1) + B_Tensor(2,2)
            ! dev_B(:,:) = B_Tensor(:,:) - trace_B*I_Tensor(:,:)/2.0d0
            ! Te_Tensor(:,:) = GeN*dev_B(:,:)

            P_Tensor(:,:) = Tve_Tensor(:,:) + SvN*G_dot(:,:) - P*I_Tensor(:,:) &
                + Te_Tensor(:,:) + Ta_Tensor(:,:)

            div_P(:) = dPdXi(:)
            Mom(:) = div_P(:)

            Continuity = traceTensor(GV)

            !Stabilization
            Ha = GveN**2+doubleDotProduct(Tve_Tensor,Tve_Tensor) !original
            if (GveN < 1.0d-8) then
                Ha = GeN**2+doubleDotProduct(Te_Tensor,Te_Tensor)
            end if
            Ha = Ha/(UlN**2+doubleDotProduct(G_dot,G_dot)) !original
            ! Ha = Ha/(EdN**2+doubleDotProduct(G_dot,G_dot))
            Ha = sqrt(Ha)

            tpspg = (SvN/Helem**2)**2 + (Ha/Helem**2)**2 !original
            tpspg = 1.0d0/sqrt(tpspg)

            tlsic = 1.0d0/UlN !original
            ! tlsic = 1.0d0/EdN
            ! tlsic = (Helem**2)/tpspg
            ! tlsic = GveN
            ! tlsic = GeN

            LOOP_RESIDUALS_f:do inod = 1, NBF_3d

                BIFN = bfn(inod)
                DBIX = dfdxi(inod,1)
                DBIY = dfdxi(inod,2)
                DBIZ = dfdxi(inod,3)

                TERM_RES(1) = P_Tensor(1,1)*DBIX + P_Tensor(2,1)*DBIY + P_Tensor(3,1)*DBIZ &
                            + 0.0d0!tlsic*(Continuity)*DBIX

                TERM_RES(2) = P_Tensor(1,2)*DBIX + P_Tensor(2,2)*DBIY + P_Tensor(3,2)*DBIZ &
                            + 0.0d0!tlsic*(Continuity)*DBIY

                TERM_RES(3) = P_Tensor(1,3)*DBIX + P_Tensor(2,3)*DBIY + P_Tensor(3,3)*DBIZ &
                            + 0.0d0!tlsic*(Continuity)*DBIZ

                TERM_RES(4) = (Continuity)*BIFN &
                            + tpspg*(Mom(1)*DBIX+Mom(2)*DBIY+Mom(3)*DBIZ) !+ (detS-1.0d0)*BIFN

                TEMP_RES(inod,:) = TEMP_RES(inod,:) + TERM_RES(:)*WET 

            end do LOOP_RESIDUALS_f

        end do LOOP_GAUSS

        if ( STORE ) then
            NM = NM_f(iel,1:NBF_3d)
            call storeResidual( TEMP_RES, NM, NBF_3d, NEQ_f, B_f, NUNKNOWNS_f )
        end if

    end subroutine bulkResidual

    !---------------------------------------------------------------------

    subroutine bulkJacobian( iel, TEMP_TL, TEMP_TD, TEMP_RES )

        use ELEMENTS_MODULE, only: NBF_3d, NEQ_f, NEQ_s, NUNKNOWNS_f
        use ENUMERATION_MODULE, only: NM_f
        use CSR_STORAGE, only: A_f, IA_f, CSR_f, NZ_f
        use Tools_mod, only: perturbVariable, storeJacobian

        implicit none

        integer, intent(in) :: iel
        real(8), dimension(NBF_3d,NEQ_f), intent(inout) :: TEMP_TL
        real(8), dimension(NEQ_s) :: TEMP_TD
        real(8), dimension(NBF_3d,NEQ_f), intent(in) :: TEMP_RES

        integer :: inod, jnod, ieq
        real(8) :: EPS_JAC

        integer, dimension(NBF_3d) :: NM
        integer, dimension(NBF_3d*NBF_3d) :: CSC

        real(8), dimension(NBF_3d,NEQ_f) :: dTEMP_RES

        real(8), dimension(NBF_3d,NBF_3d, NEQ_f, NEQ_f) :: TEMP_JAC

        TEMP_JAC  = 0.0d0

        do jnod = 1, NBF_3d
            do ieq = 1, NEQ_f
                EPS_JAC = perturbVariable( TEMP_TL(jnod,ieq) )
                TEMP_TL(jnod,ieq) = TEMP_TL(jnod,ieq) + EPS_JAC
                call bulkResidual( iel, TEMP_TL, TEMP_TD, dTEMP_RES, .false. )  
                TEMP_TL(jnod,ieq) = TEMP_TL(jnod,ieq) - EPS_JAC
                do inod = 1, NBF_3d
                    TEMP_JAC(inod,jnod,ieq,:) = ( dTEMP_RES(inod,:) - TEMP_RES(inod,:) )/EPS_JAC
                end do
            end do
        end do

        NM  = NM_f(iel,1:NBF_3d)
        CSC = CSR_f(iel,1:NBF_3d*NBF_3d)

        call storeJacobian&
            (TEMP_JAC, NBF_3d, NBF_3d, NEQ_f, NEQ_f, NM, IA_f, NUNKNOWNS_f+1,&
            CSC, NBF_3d*NBF_3d, A_f, NZ_f)

    end subroutine bulkJacobian

    !---------------------------------------------------------------------

    subroutine STRESS_STORAGE( iel )

        use ELEMENTS_MODULE, only: NBF_3d, NEQ_f, NEQ_s
        use ENUMERATION_MODULE, only: NM_MESH
        use GLOBAL_ARRAYS_MODULE, only: TL, TD
        use Tools_mod, only: determinant3D
        use Stress_mod, only: DOMI_STRESS

        implicit none

        integer, intent(in)  :: iel

        real(8), dimension(3,3) :: C_Tensor, S_Tensor, A_Tensor
        real(8) :: detC, detS, detA

        real(8), dimension(NBF_3d,NEQ_f) :: TEMP_TL
        real(8), dimension(NEQ_s) :: TEMP_TD

        TEMP_TL(:,:) = TL(NM_MESH(iel,:),:)
        TEMP_TD(:) = TD(iel,:)

        call DOMI_STRESS( iel, TEMP_TL, TEMP_TD )

        TD(iel,:) = TEMP_TD(:)

        !C Tensor
        C_Tensor(1,1) = TEMP_TD(1) ; C_Tensor(1,2) = TEMP_TD(2) ; C_Tensor(1,3) = TEMP_TD(3)
        C_Tensor(2,1) = TEMP_TD(2) ; C_Tensor(2,2) = TEMP_TD(4) ; C_Tensor(2,3) = TEMP_TD(5)
        C_Tensor(3,1) = TEMP_TD(3) ; C_Tensor(3,2) = TEMP_TD(5) ; C_Tensor(3,3) = TEMP_TD(6)

        detC = determinant3D(C_Tensor)
        if (detC < 0.0d0) write(*,'(a,i8,es12.4)') 'detC < 0 at iel:', iel, detC

        !S Tensor
        S_Tensor(1,1) = TEMP_TD(7) ; S_Tensor(1,2) = TEMP_TD(8) ; S_Tensor(1,3) = TEMP_TD(9)
        S_Tensor(2,1) = TEMP_TD(8) ; S_Tensor(2,2) = TEMP_TD(10) ; S_Tensor(2,3) = TEMP_TD(11)
        S_Tensor(3,1) = TEMP_TD(9) ; S_Tensor(3,2) = TEMP_TD(11) ; S_Tensor(3,3) = TEMP_TD(12)

        detS = determinant3D(S_Tensor)
        ! if ( abs(detS-1.0d0) > 1.0d-8) write(*,'(a,i8,es12.4)') 'detS /= 1 at iel:', iel, detS
        if (detS < 0.0d0) write(*,'(a,i8,es12.4)') 'detS < 0 at iel:', iel, detS

        !A Tensor
        A_Tensor(1,1) = TEMP_TD(13) ; A_Tensor(1,2) = TEMP_TD(14) ; A_Tensor(1,3) = TEMP_TD(15)
        A_Tensor(2,1) = TEMP_TD(14) ; A_Tensor(2,2) = TEMP_TD(16) ; A_Tensor(2,3) = TEMP_TD(17)
        A_Tensor(3,1) = TEMP_TD(15) ; A_Tensor(3,2) = TEMP_TD(17) ; A_Tensor(3,3) = TEMP_TD(18)

        detA = determinant3D(A_Tensor)
        ! if ( abs(detA-1.0d0) > 1.0d-8) write(*,'(a,i8,es12.4)') 'detA /= 1 at iel:', iel, detA
        if (detA < 0.0d0) write(*,'(a,i8,es12.4)') 'detA < 0 at iel:', iel, detA

    end subroutine STRESS_STORAGE

    !----------------------------------------------------------------------

    subroutine DIRICHLET_BOUNDARY_CONDITIONS( NOPT )

        use TIME_INTEGRATION_MODULE, only: Total_time
        use PHYSICAL_MODULE, only: X0, Y0, Z0, STRAIN
        use ELEMENTS_MODULE, only: NODTOL
        use ENUMERATION_MODULE, only: GNTR
        use FLOW_ARRAYS_MODULE, only: Bc_f, Bce_f, B_f
        use GLOBAL_ARRAYS_MODULE, only: TL
        use MESH_MODULE, only: Xm, Ym, Zm, EPS_MESH
        use CSR_STORAGE, only: A_f, IA_f, DA_f
        use BOUNDARY_ENUMERATION_MODULE, only: BND_Symmetry, Symmetry

        implicit none

        integer, intent(in) :: NOPT

        integer :: irow, iL, iU, jcol, inod, ibnd, i
                
        select case(NOPT)

        case(1)

            Bc_f  = .false.
            Bce_f = .false.

            do inod = 1, NODTOL

                IROW = GNTR(inod)

                do i = 1, size(BND_Symmetry)
                    ibnd = BND_Symmetry(i)

                    select case (Symmetry(i))

                    case ('X')

                        if ( abs(Xm(inod)-0.0d0)<EPS_MESH ) then
                            B_f(IROW) = TL(inod,1) - Xm(inod)
                            Bc_f(IROW) = .true.
                        end if

                    case ('Y')

                        if ( abs(Ym(inod)-0.0d0)<EPS_MESH ) then
                            B_f(IROW+1) = TL(inod,2) - Ym(inod)
                            Bc_f(IROW+1) = .true.
                        end if

                    case ('Z')

                        if ( abs(Zm(inod)-0.0d0)<EPS_MESH ) then
                            B_f(IROW+2) = TL(inod,3) - Zm(inod)
                            Bc_f(IROW+2) = .true.
                        end if

                    case default
                        write(*,'(a)') 'Wrong Symmetry condition in DIRICHLET_BOUNDARY_CONDITIONS!'
                        stop

                    end select

                end do

                if ( abs(Zm(inod)-0.0d0)<EPS_MESH ) then

                    B_f(IROW) = TL(inod,1) - Xm(inod)
                    Bc_f(IROW) = .true.

                    B_f(IROW+1) = TL(inod,2) - Ym(inod)
                    Bc_f(IROW+1) = .true.

                    B_f(IROW+2) = TL(inod,3) - Zm(inod)
                    Bc_f(IROW+2) = .true.

                end if

                if ( abs(Zm(inod)-Z0)<EPS_MESH ) then

                    B_f(IROW) = TL(inod,1) - Xm(inod)
                    Bc_f(IROW) = .true.

                    B_f(IROW+1)  = TL(inod,2) - Ym(inod)
                    Bc_f(IROW+1) = .true.

                    B_f(IROW+2) = TL(inod,3) - STRAIN(Total_time)
                    Bc_f(IROW+2) = .true.

                end if

            end do

        case(2)

            do IROW = 1, size(Bc_f)
                if (Bc_f(IROW)) then
                    iL = IA_f(IROW)
                    iU = IA_f(IROW+1)-1
                    A_f(iL:iU) = 0.0d0

                    jcol = DA_f(IROW)
                    A_f(jcol)  = 1.0d0   
                end if

                if (Bce_f(IROW)) then
                    iL = IA_f(IROW)
                    iU = IA_f(IROW+1)-1

                    jcol = DA_f(IROW)+1
                    A_f(jcol)  = -1.0d0!/tan(SaN)
                end if
            end do

        case default

            write(*,*)'Incorrect option in Dirichlet BCs: ', NOPT
            stop

        end select

    end subroutine DIRICHLET_BOUNDARY_CONDITIONS

    !---------------------------------------------------------------------

    subroutine INTERFACE_FORCE_RESIDUAL_f( iel, ifc, TEMP_TL, TEMP_RES, STORE )

        use TIME_INTEGRATION_MODULE, only: INCREMENT, dt
        use PHYSICAL_MODULE, only: PaN, SvN
        use ELEMENTS_MODULE, only: NBF_3d, NEL_3d, NEQ_f, NUNKNOWNS_f, NEQ_s
        use GAUSS_MODULE, only: WO_2d, NGAUSS_2d, BFN_F, DFDC_F, DFDE_F, DFDS_F
        use ENUMERATION_MODULE, only: NM_MESH, NM_f
        use GLOBAL_ARRAYS_MODULE, only: TLo, TLb, TD
        use FLOW_ARRAYS_MODULE, only: B_f
        use MESH_MODULE, only: Xm, Ym, Zm
        use Stress_mod, only: DOMI_STRESS, VE_tensor, Elastic_tensor, Orientation_tensor
        use Tools_mod, only: storeResidual, basis3D

        implicit none

        integer, intent(in)  :: iel, ifc
        real(8), dimension(NBF_3d, NEQ_f), intent(in)  :: TEMP_TL
        real(8), dimension(NBF_3d, NEQ_f), intent(out) :: TEMP_RES
        logical, intent(in)  :: STORE

        integer, parameter, dimension(3,3) :: I_Tensor = [1,0,0,0,1,0,0,0,1]
        integer :: ig, inod, gnod, i
        real(8) :: WET, P, CJAC, dS

        real(8) :: bifn
        real(8), dimension(3) :: normal
        real(8), dimension(3,3) :: C_Tensor, Tve_Tensor, G_dot, GV, GVT, P_Tensor, A_Tensor, Ta_Tensor
        real(8), dimension(3,3) :: GU, GUT, F_Tensor, FT_Tensor, B_Tensor, Te_Tensor, dev_B, S_Tensor

        integer, dimension(NBF_3d) :: NM
        real(8), dimension(NEQ_f)  :: TERM_RES
        real(8), dimension(NEQ_s)  :: TEMP_TD

        real(8), dimension(3) :: Xi_gs, dt_term
        real(8), dimension(NBF_3d) :: bfn
        real(8), dimension(3,3) :: dXidci
        real(8), dimension(NBF_3d,3) :: Xi_el, dfdci, dfdxi

        Xi_el(:,:) = TEMP_TL(:,1:3)
        ! Xi_el(:,1) = Xm(NM_MESH(iel,:))
        ! Xi_el(:,2) = Ym(NM_MESH(iel,:))
        ! Xi_el(:,3) = Zm(NM_MESH(iel,:))

        TEMP_TD = TD(iel,:)

        call DOMI_STRESS( iel, TEMP_TL,  TEMP_TD )

        C_Tensor(1,1) = TEMP_TD(1) ; C_Tensor(1,2) = TEMP_TD(2) ; C_Tensor(1,3) = TEMP_TD(3)
        C_Tensor(2,1) = TEMP_TD(2) ; C_Tensor(2,2) = TEMP_TD(4) ; C_Tensor(2,3) = TEMP_TD(5)
        C_Tensor(3,1) = TEMP_TD(3) ; C_Tensor(3,2) = TEMP_TD(5) ; C_Tensor(3,3) = TEMP_TD(6)

        S_Tensor(1,1) = TEMP_TD(7) ; S_Tensor(1,2) = TEMP_TD(8)  ; S_Tensor(1,3) = TEMP_TD(9)
        S_Tensor(2,1) = TEMP_TD(8) ; S_Tensor(2,2) = TEMP_TD(10) ; S_Tensor(2,3) = TEMP_TD(11)
        S_Tensor(3,1) = TEMP_TD(9) ; S_Tensor(3,2) = TEMP_TD(11) ; S_Tensor(3,3) = TEMP_TD(12)

        A_Tensor(1,1) = TEMP_TD(13) ; A_Tensor(1,2) = TEMP_TD(14) ; A_Tensor(1,3) = TEMP_TD(15)
        A_Tensor(2,1) = TEMP_TD(14) ; A_Tensor(2,2) = TEMP_TD(16) ; A_Tensor(2,3) = TEMP_TD(17)
        A_Tensor(3,1) = TEMP_TD(15) ; A_Tensor(3,2) = TEMP_TD(17) ; A_Tensor(3,3) = TEMP_TD(18)

        call VE_tensor(C_Tensor, Tve_Tensor)
        call Elastic_tensor(S_Tensor, Te_Tensor)
        call Orientation_tensor(A_Tensor, Ta_Tensor)

        TEMP_RES = 0.0d0

        LOOP_GAUSS: do ig = 1, NGAUSS_2d

            bfn(:) = BFN_F(:,ig,ifc)
            dfdci(:,1) = DFDC_F(:,ig,ifc)
            dfdci(:,2) = DFDE_F(:,ig,ifc)
            dfdci(:,3) = DFDS_F(:,ig,ifc)

            call basis3D(Xi_el, bfn, dfdci, dfdxi, Xi_gs, dXidci, CJAC)

            associate(dXdC => dXidci(1,1), dXdE => dXidci(1,2), dXdS => dXidci(1,3), &
                dYdC => dXidci(2,1), dYdE => dXidci(2,2), dYdS => dXidci(2,3), &
                dZdC => dXidci(3,1), dZdE => dXidci(3,2), dZdS => dXidci(3,3))

                select case(ifc)
                case(1)
                    normal(1) =  (dYdE*dZdS-dZdE*dYdS)
                    normal(2) = -(dXdE*dZdS-dZdE*dXdS)
                    normal(3) =  (dXdE*dYdS-dYdE*dXdS)

                case(2)
                    normal(1) = -(dYdC*dZdS-dZdC*dYdS)
                    normal(2) =  (dXdC*dZdS-dZdC*dXdS)
                    normal(3) = -(dXdC*dYdS-dYdC*dXdS)

                case(3)
                    normal(1) =  (dYdC*dZdE-dZdC*dYdE)
                    normal(2) = -(dXdC*dZdE-dZdC*dXdE)
                    normal(3) =  (dXdC*dYdE-dYdC*dXdE)

                case(4)
                    normal(1) = -((dYdC-dYdS)*(dZdE-dZdS)-(dZdC-dZdS)*(dYdE-dYdS))
                    normal(2) =  ((dXdC-dXdS)*(dZdE-dZdS)-(dZdC-dZdS)*(dXdE-dXdS))
                    normal(3) = -((dXdC-dXdS)*(dYdE-dYdS)-(dYdC-dYdS)*(dXdE-dXdS))

                end select

            end associate

            dS = sqrt(normal(1)**2+normal(2)**2+normal(3)**2)

            WET = WO_2d(ig)*dS

            normal(:) = normal(:)/dS

            P = 0.0d0
            GV = 0.0d0
            do inod = 1, NBF_3d

                gnod = NM_MESH(iel,inod)

                P = P + TEMP_TL(inod,4)*bfn(inod)

                if (Increment == 1) then
                    dt_term(:) = (TEMP_TL(inod,1:3)-TLo(gnod,1:3))/dt
                else
                    dt_term(:) = (3.0d0*TEMP_TL(inod,1:3)-4.0d0*TLo(gnod,1:3)+TLb(gnod,1:3))/(2.0d0*dt)
                end if
                do i = 1, 3
                    GV(:,i) = GV(:,i) + dt_term(i)*dfdxi(inod,:)
                end do

            end do

            GVT(:,:) = transpose(GV)
            G_dot(:,:) = GV(:,:) + GVT(:,:)

            P_Tensor(:,:) = Tve_Tensor(:,:) + SvN*G_dot(:,:) - P*I_Tensor(:,:) &
                + Te_Tensor(:,:) + Ta_Tensor(:,:)

            LOOP_RESIDUALS_f:DO inod = 1, NBF_3d

                bifn = bfn(inod)

                TERM_RES(1) = ( dot_product(normal,P_Tensor(1,:))+PaN )*bifn
                TERM_RES(2) = ( dot_product(normal,P_Tensor(2,:))+PaN )*bifn
                TERM_RES(3) = ( dot_product(normal,P_Tensor(3,:))+PaN )*bifn

                TERM_RES(4) = 0.0d0

                TEMP_RES(inod,:) = TEMP_RES(inod,:) + TERM_RES(:)*WET 

            end do LOOP_RESIDUALS_f

        end do LOOP_GAUSS

        if ( STORE ) then
            NM = NM_f(iel,1:NBF_3d)
            call storeResidual( TEMP_RES, NM, NBF_3d, NEQ_f, B_f, NUNKNOWNS_f )
        end if

    end subroutine INTERFACE_FORCE_RESIDUAL_f

    !---------------------------------------------------------------------

    subroutine INTERFACE_FORCE_JACOBIAN_f( iel, ifc, TEMP_TL, TEMP_RES )

        use ELEMENTS_MODULE, only: NBF_3d, NEQ_f, NUNKNOWNS_f
        use ENUMERATION_MODULE, only: NM_f
        use CSR_STORAGE, only: A_f, IA_f, CSR_f, NZ_f
        use Tools_mod, only: perturbVariable, storeJacobian

        implicit none

        integer, intent(in) :: iel, ifc
        real(8), dimension(NBF_3d,NEQ_f), intent(inout) :: TEMP_TL
        real(8), dimension(NBF_3d,NEQ_f), intent(in)    :: TEMP_RES

        integer :: II, JJ, IW, JW, I, J, INOD, IEQ, JNOD, JEQ
        integer :: IROW, JCOL, ICOL, IAD, L
        real(8) :: EPS_JAC

        integer, dimension(NBF_3d)        :: NM
        integer, dimension(NBF_3d*NBF_3d) :: CSC

        real(8), dimension(NBF_3d,NEQ_f)  :: dTEMP_RES

        real(8), dimension(NBF_3d,NBF_3d, NEQ_f, NEQ_f) :: TEMP_JAC

        TEMP_JAC  = 0.0d0

        do JW = 1, NBF_3d
            do JEQ = 1, NEQ_f
                EPS_JAC = perturbVariable( TEMP_TL(JW,JEQ) )    
                TEMP_TL(JW,JEQ) = TEMP_TL(JW,JEQ) + EPS_JAC
                call INTERFACE_FORCE_RESIDUAL_f( iel, ifc, TEMP_TL, dTEMP_RES, .false. )  
                TEMP_TL(JW,JEQ) = TEMP_TL(JW,JEQ) - EPS_JAC
                do IW = 1, NBF_3d
                    TEMP_JAC(IW,JW,JEQ,1:NEQ_f) = ( dTEMP_RES(IW,1:NEQ_f) - TEMP_RES(IW,1:NEQ_f) )/EPS_JAC
                end do
            end do
        end do

        NM  = NM_f(iel,1:NBF_3d)
        CSC = CSR_f(iel,1:NBF_3d*NBF_3d)

        call storeJacobian&
            (TEMP_JAC, NBF_3d, NBF_3d, NEQ_f, NEQ_f, NM, IA_f, NUNKNOWNS_f+1,&
            CSC, NBF_3d*NBF_3d, A_f, NZ_f)

    end subroutine INTERFACE_FORCE_JACOBIAN_f

end module Equations_mod
