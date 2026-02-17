
module PostProcess_mod

    implicit none

    contains

    subroutine postProcess(Increment)

        use PHYSICAL_MODULE, only: Param
        use TIME_INTEGRATION_MODULE, only: Total_time, Initial_time
        use GLOBAL_ARRAYS_MODULE, only: TL, TLo, TLb, TD, TDo, TDb
        use Tools_mod, only: plotGraphs

        implicit none

        integer, intent(in) :: Increment

        real(8), parameter :: eps = 1.0d-8, Dt_write = 20.0d0!0.01d0
        real(8), save :: twrite0 = Initial_time
        real(8) :: e_eng, stress_eng
        logical :: check
        character(100) :: fn, str, fn_temp
        character(len=:), allocatable :: str1

        if (Increment == 1) twrite0 = Initial_time

        check = (Total_time - twrite0 > Dt_write - eps)
        ! if (Total_time > 10.0d0 + eps) check = (Total_time - twrite0 > 10.0d0*Dt_write - eps)
        ! if (Total_time > 1000.0d0 + eps) check = (Total_time - twrite0 > 50.0d0*Dt_write - eps)
        ! if (Total_time > 5000.0d0 + eps) check = (Total_time - twrite0 > 100.0d0*Dt_write - eps)
        ! if (Total_time > 10000.0d0 + eps) check = (Total_time - twrite0 > 500.0d0*Dt_write - eps)
        ! ! check = ( mod(Increment,100) == 0 )

        check = check .or. (Increment == 1)

        if (check) then

            twrite0 = Total_time

            if (Increment == 1) twrite0 = 0.0d0

            call writeTecplotFile(Total_time)

            e_eng = specimenSize()

            stress_eng = stress_strain(e_eng, check)

            ! check = plotGraphs(Increment)

            fn_temp = 'out/SOL/SOL_'
            if (associated(Param%p)) then
                write(str,'(f12.4)') Param%p
                str1 = trim(adjustl(str))
                fn_temp = trim(adjustl(fn_temp))//trim(adjustl(Param%name))//'_'//str1
                deallocate(str1)
            end if

            write(str,'(f12.4)') Total_time
            str1 = trim(adjustl(str))

            fn = trim(adjustl(fn_temp))//'_'//str1
            call writeSolution(fn, TL, TD)

            fn = trim(adjustl(fn_temp))//'_'//str1//'_o'
            call writeSolution(fn, TLo, TDo)

            fn = trim(adjustl(fn_temp))//'_'//str1//'_b'
            call writeSolution(fn, TLb, TDb)

            deallocate(str1)

        end if

    end subroutine postProcess

    !---------------------------------------------------------------------

    function specimenSize() result(e_eng)

        use TIME_INTEGRATION_MODULE, only: Total_time
        use PHYSICAL_MODULE, only: Z0, Param
        use ELEMENTS_MODULE, only: NEL_3d
        use BOUNDARY_ENUMERATION_MODULE, only: NBE

        implicit none

        real(8) :: e_eng

        integer :: iel, ibnd
        real(8) :: Xmin, z_Xmin, Ymin, z_Ymin, zmax, e_true
        character(100) :: fn, str
        character(len=:), allocatable :: str1

        Xmin = huge(0.0d0) ; Ymin = huge(0.0d0) ; zmax = 0.0d0
        do iel = 1, NEL_3d

            !X=X0
            ibnd = 2
            if (NBE(iel)%KOE(ibnd)) then
                call filamentXsize( iel, NBE(iel)%nodes(:,ibnd), Xmin, z_Xmin, zmax )
            end if

            !Y=Yi
            ibnd = 4
            if (NBE(iel)%KOE(ibnd)) then
                call filamentYsize( iel, NBE(iel)%nodes(:,ibnd), Ymin, z_Ymin, zmax )
            end if

        end do

        e_eng = (zmax-Z0)/Z0
        ! e_true = log(zmax/Z0)
        ! e_true = log(1.0d0+e_eng)

        fn = 'out/Specimen_size'
        if (associated(Param%p)) then
            write(str,'(f12.4)') Param%p
            str1 = trim(adjustl(str))
            fn = trim(adjustl(fn))//'_'//trim(adjustl(Param%name))//'_'//str1
            deallocate(str1)
        end if
        fn = trim(adjustl(fn))//'.DAT'

        open(20,file=fn,position='append')
        write(20,'(*(f16.6))') Total_time, Xmin, z_Xmin, Ymin, z_Ymin, zmax
        close(20)

    end function specimenSize

    !---------------------------------------------------------------------

    subroutine filamentXsize( iel, nodes, Xmin, z_Xmin, zmax )

        use GLOBAL_ARRAYS_MODULE, only: TL

        implicit none

        integer, intent(in) :: iel
        integer, dimension(:), intent(in) :: nodes
        real(8), intent(inout) :: Xmin, z_Xmin, zmax

        integer :: inod, gnod
        real(8) :: x, z

        do inod = 1, size(nodes)

            gnod = nodes(inod)

            x = TL(gnod,1)
            z = TL(gnod,3)

            if (x < Xmin) then
                Xmin = x
                z_Xmin = z
            end if

            if (z > zmax) zmax = z

        end do

    end subroutine filamentXsize

    !---------------------------------------------------------------------

    subroutine filamentYsize( iel, nodes, Ymin, z_Ymin, zmax )

        use GLOBAL_ARRAYS_MODULE, only: TL

        implicit none

        integer, intent(in) :: iel
        integer, dimension(:), intent(in) :: nodes
        real(8), intent(inout) :: Ymin, z_Ymin, zmax

        integer :: inod, gnod
        real(8) :: y, z

        do inod = 1, size(nodes)

            gnod = nodes(inod)

            y = TL(gnod,2)
            z = TL(gnod,3)

            if (y < Ymin) then
                Ymin = y
                z_Ymin = z
            end if

            if (z > zmax) zmax = z

        end do

    end subroutine filamentYsize
    
    !---------------------------------------------------------------------

    function stress_strain(e_eng, output) result(stress_eng)

        use TIME_INTEGRATION_MODULE, only: Total_time
        use PHYSICAL_MODULE, only: X0, Y0, Z0, Param
        use ELEMENTS_MODULE, only: NEL_3d
        use BOUNDARY_ENUMERATION_MODULE, only: NBE, Nbd_Symmetry

        implicit none

        real(8), intent(in) :: e_eng
        logical, intent(in) :: output
        real(8) :: stress_eng

        integer :: iel, ibnd
        real(8) :: force, force_el, e_true, stress_true
        character(100) :: fn, str
        character(len=:), allocatable :: str1

        force = 0.0d0
        do iel = 1, NEL_3d

            !Z=Z0
            ibnd = 6
            if (NBE(iel)%KOE(ibnd)) then
                force_el = force2D( iel, NBE(iel)%BFE(ibnd) )
                force = force + 2.0d0*Nbd_Symmetry*force_el !in N
            end if

        end do

        ! e_true = log(zmax/Z0)
        e_true = log(1.0d0+e_eng)

        stress_eng = force/(2.0d0*Nbd_Symmetry*X0*Y0) !in MPa
        stress_true = stress_eng*(1.0d0+e_eng)

        if (output) then

            fn = 'out/Stress'
            if (associated(Param%p)) then
                write(str,'(f12.4)') Param%p
                str1 = trim(adjustl(str))
                fn = trim(adjustl(fn))//'_'//trim(adjustl(Param%name))//'_'//str1
                deallocate(str1)
            end if
            fn = trim(adjustl(fn))//'.DAT'

            open(21,file=fn,position='append')
            write(21,'(*(es20.8))') Total_time, e_eng*100.0d0, e_true*100.0d0, stress_eng, stress_true
            close(21)
        end if

    end function stress_strain
    
    !---------------------------------------------------------------------

    function stress_strain_true(e_true, output) result(stress_true)

        use TIME_INTEGRATION_MODULE, only: Total_time
        use PHYSICAL_MODULE, only: X0, Y0, Z0, Param
        use ELEMENTS_MODULE, only: NEL_3d
        use BOUNDARY_ENUMERATION_MODULE, only: NBE, Nbd_Symmetry

        implicit none

        real(8), intent(in) :: e_true
        logical, intent(in) :: output
        real(8) :: stress_true

        integer :: iel, ibnd
        real(8) :: force, force_el, e_eng, stress_eng
        character(100) :: fn, str
        character(len=:), allocatable :: str1

        force = 0.0d0
        do iel = 1, NEL_3d

            !Z=Z0
            ibnd = 6
            if (NBE(iel)%KOE(ibnd)) then
                force_el = force2D( iel, NBE(iel)%BFE(ibnd) )
                force = force + 2.0d0*Nbd_Symmetry*force_el !in N
            end if

        end do

        e_eng = exp(e_true)-1.0d0

        stress_eng = force/(2.0d0*Nbd_Symmetry*X0*Y0) !in MPa
        stress_true = stress_eng*(1.0d0+e_eng)

        if (output) then

            fn = 'out/Stress'
            if (associated(Param%p)) then
                write(str,'(f12.4)') Param%p
                str1 = trim(adjustl(str))
                fn = trim(adjustl(fn))//'_'//trim(adjustl(Param%name))//'_'//str1
                deallocate(str1)
            end if
            fn = trim(adjustl(fn))//'.DAT'

            open(21,file=fn,position='append')
            write(21,'(*(es20.8))') Total_time, e_eng*100.0d0, e_true*100.0d0, stress_eng, stress_true
            close(21)
        end if

    end function stress_strain_true

    !---------------------------------------------------------------------

    function force2D(iel, ifc) result(force_el)

        use TIME_INTEGRATION_MODULE, only: Increment, dt
        use PHYSICAL_MODULE, only: SvN
        use ELEMENTS_MODULE, only: NBF_3d, NEQ_f
        use ENUMERATION_MODULE, only: NM_MESH
        use GLOBAL_ARRAYS_MODULE, only: TL, TLo, TLb, TD
        ! use MESH_MODULE, only: Xm, Ym, Zm
        use GAUSS_MODULE, only: WO_2d, NGAUSS_2d, BFN_F, DFDC_F, DFDE_F, DFDS_F
        use Stress_mod, only: VE_tensor, Elastic_tensor, Orientation_tensor
        use Tools_mod, only: basis3D

        implicit none

        integer, intent(in) :: iel, ifc
        real(8) :: force_el

        integer :: ig, inod, gnod, i
        real(8) :: WET, P, CJAC, dS, nx, ny, nz
        integer, parameter, dimension(3,3) :: I_Tensor = [1,0,0,0,1,0,0,0,1]
        real(8), dimension(3,3) :: C_Tensor, Tve_Tensor, GV, GVT, G_dot, S_Tensor, &
            Te_Tensor, P_Tensor, A_Tensor, Ta_Tensor
        real(8), dimension(NBF_3d, NEQ_f) :: TEMP_TL
        real(8), dimension(3) :: Xi_gs, dt_term
        real(8), dimension(NBF_3d) :: bfn
        real(8), dimension(3,3) :: dXidci
        real(8), dimension(NBF_3d,3) :: Xi_el, dfdci, dfdxi

        TEMP_TL = TL(NM_MESH(iel,:),:)

        Xi_el(:,:) = TEMP_TL(:,1:3)
        ! Xi_el(:,1) = Xm(NM_MESH(iel,:))
        ! Xi_el(:,2) = Ym(NM_MESH(iel,:))
        ! Xi_el(:,3) = Zm(NM_MESH(iel,:))

        C_Tensor(1,1) = TD(iel,1) ; C_Tensor(1,2) = TD(iel,2) ; C_Tensor(1,3) = TD(iel,3)
        C_Tensor(2,1) = TD(iel,2) ; C_Tensor(2,2) = TD(iel,4) ; C_Tensor(2,3) = TD(iel,5)
        C_Tensor(3,1) = TD(iel,3) ; C_Tensor(3,2) = TD(iel,5) ; C_Tensor(3,3) = TD(iel,6)
        call VE_tensor(C_Tensor,Tve_Tensor)

        S_Tensor(1,1) = TD(iel,7) ; S_Tensor(1,2) = TD(iel,8)  ; S_Tensor(1,3) = TD(iel,9)
        S_Tensor(2,1) = TD(iel,8) ; S_Tensor(2,2) = TD(iel,10) ; S_Tensor(2,3) = TD(iel,11)
        S_Tensor(3,1) = TD(iel,9) ; S_Tensor(3,2) = TD(iel,11) ; S_Tensor(3,3) = TD(iel,12)
        call Elastic_tensor(S_Tensor, Te_Tensor)

        A_Tensor(1,1) = TD(iel,13) ; A_Tensor(1,2) = TD(iel,14) ; A_Tensor(1,3) = TD(iel,15)
        A_Tensor(2,1) = TD(iel,14) ; A_Tensor(2,2) = TD(iel,16) ; A_Tensor(2,3) = TD(iel,17)
        A_Tensor(3,1) = TD(iel,15) ; A_Tensor(3,2) = TD(iel,17) ; A_Tensor(3,3) = TD(iel,18)
        call Orientation_tensor(A_Tensor, Ta_Tensor)

        force_el = 0.0d0
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
                    nx =   (dYdE*dZdS-dZdE*dYdS)
                    ny = - (dXdE*dZdS-dZdE*dXdS)
                    nz =   (dXdE*dYdS-dYdE*dXdS)

                case(2)
                    nx = - (dYdC*dZdS-dZdC*dYdS)
                    ny =   (dXdC*dZdS-dZdC*dXdS)
                    nz = - (dXdC*dYdS-dYdC*dXdS)

                case(3)
                    nx =   (dYdC*dZdE-dZdC*dYdE)
                    ny = - (dXdC*dZdE-dZdC*dXdE)
                    nz =   (dXdC*dYdE-dYdC*dXdE)

                case(4)
                    nx = - ((dYdC-dYdS)*(dZdE-dZdS)-(dZdC-dZdS)*(dYdE-dYdS))
                    ny =   ((dXdC-dXdS)*(dZdE-dZdS)-(dZdC-dZdS)*(dXdE-dXdS))
                    nz = - ((dXdC-dXdS)*(dYdE-dYdS)-(dYdC-dYdS)*(dXdE-dXdS))

                end select

            end associate

            dS = sqrt(nx**2+ny**2+nz**2)

            WET = WO_2d(ig)*dS

            nx = nx/dS
            ny = ny/dS
            nz = nz/dS

            GV(:,:) = 0.0d0
            P = 0.0d0
            do inod = 1, NBF_3d

                gnod = NM_MESH(iel,inod)

                P = P + TEMP_TL(inod,4)*bfn(inod)

                dt_term(:) = (3.0d0*TEMP_TL(inod,1:3)-4.0d0*TLo(gnod,1:3)+TLb(gnod,1:3))/(2.0d0*dt)
                if (Increment == 1) then
                    dt_term(:) = (TEMP_TL(inod,1:3)-TLo(gnod,1:3))/dt
                end if
                do i = 1, 3
                    GV(:,i) = GV(:,i) + dt_term(i)*dfdxi(inod,:)
                end do

            end do

            GVT(:,:) = transpose(GV)
            G_dot(:,:) = GV(:,:) + GVT(:,:)

            P_Tensor(:,:) = Tve_Tensor(:,:) + SvN*G_dot(:,:) - P*I_Tensor(:,:) &
                + Te_Tensor(:,:) + Ta_Tensor(:,:)

            force_el = force_el + P_Tensor(3,3)*WET

        end do LOOP_GAUSS

    end function force2D

    !---------------------------------------------------------------------

    function stressProjection() result(TDg)

        use ELEMENTS_MODULE, only: NBF_3d, NODTOL, NEQ_s
        use ENUMERATION_MODULE, only: NM_MESH, EsuN1, EsuN2
        use GLOBAL_ARRAYS_MODULE, only: TD
        use MESH_MODULE, only: Xm, Ym, Zm, EPS_MESH

        implicit none

        real(8), dimension(NODTOL,NEQ_s) :: TDg

        integer :: inod, gnod, i, j, iel
        real(8) :: Wi, SWi, Xc, Yc, Zc

        !INITIALIZE STRESS PROJECTION ARRAY
        TDg = 0.0d0
        do i = 1, NODTOL
            !ITERATE OVER  SURROUNDING ELEMENTS
            SWi = 0.0d0
            do j = EsuN2(i)+1, EsuN2(i+1)
                !# OF NEIGHBOURING ELEMENT
                iel = EsuN1(j)
                !FIND CENTROID OF K SURROUNDING ELEMENT
                Xc = 0.0d0 ; Yc = 0.0d0 ; Zc = 0.0d0
                do inod = 1, NBF_3d
                    gnod = NM_MESH(iel,inod)
                    Xc = Xc + Xm(gnod)/dble(NBF_3d)
                    Yc = Yc + Ym(gnod)/dble(NBF_3d)
                    Zc = Zc + Zm(gnod)/dble(NBF_3d)
                end do

                !DEFINE WEIGHT & UPDATE WEIGHTING SUM
                Wi = 1.0d0/sqrt((Xm(I)-Xc)**2+(Ym(I)-Yc)**2+(Zm(I)-Zc)**2)
                SWi = SWi + Wi

                !INVERSE DISTANCE INTEPOLATION
                TDg(i,:) = TDg(i,:) + Wi*TD(EsuN1(j),:)

            end do

            !DIVIDE BY SUM OF WEIGHTS
            TDg(i,:) = TDg(i,:)/SWi

            ! if ( abs(Ym(i)-0.0d0)<EPS_MESH ) then
            !     !C
            !     TDg(i,2) = 0.0d0
            !     TDg(i,5) = 0.0d0
            !     !S
            !     TDg(i,8) = 0.0d0
            !     TDg(i,11) = 0.0d0
            ! end if

        end do

    end function stressProjection

    !-----------------------------------------------------------------------

    subroutine writeTecplotFile(iter)

        use PHYSICAL_MODULE, only: Xparam, Param
        use ELEMENTS_MODULE, only: NODTOL, NEL_3d, NBF_3d
        use ENUMERATION_MODULE, only: NM_MESH
        use GLOBAL_ARRAYS_MODULE, only: TL, TDg
        use Stress_mod, only: VE_tensor, Elastic_tensor, Orientation_tensor
        use Tools_mod, only: getParameterID

        implicit none

        real(8), intent(in) :: iter

        real(8), dimension(3,3) :: S_Tensor, Te_Tensor, A_Tensor, Ta_Tensor, C_Tensor, Tve_Tensor
        real(8), dimension(NODTOL,6) :: Te_Tensor_tot, Ta_Tensor_tot, Tve_Tensor_tot
        character(len=:), allocatable :: str1
        character(100) :: str, fn
        character(500) :: variables
        integer :: ierror, i, j, ieq
        real(8) :: GveN

        GveN = Xparam(getParameterID('GveN'))

        fn = 'out/Tec/TEC_'
        if (associated(Param%p)) then
            write(str,'(f12.4)') Param%p
            str1 = trim(adjustl(str))
            fn = trim(adjustl(fn))//trim(adjustl(Param%name))//'_'//str1//'_'
            deallocate(str1)
        end if

        write(str,'(f12.4)') iter
        str1 = trim(adjustl(str))

        fn = trim(adjustl(fn))//str1//'.PLT'
        fn = trim(adjustl(fn))

        open(14, file=fn, status='unknown', action='write', iostat=ierror, position='rewind')

        if( ierror /= 0 )then
            write(*,'(a)')'Error:  PLOTS.plt'
            stop
        end if

        !PROJECT STRESSES TO CONTINUOUS ELEMENTS
        TDg(:,:) = stressProjection()

        do i = 1, NODTOL

            !C, Tve
            C_Tensor(1,1) = TDg(i,1) ; C_Tensor(1,2) = TDg(i,2) ; C_Tensor(1,3) = TDg(i,3)
            C_Tensor(2,1) = TDg(i,2) ; C_Tensor(2,2) = TDg(i,4) ; C_Tensor(2,3) = TDg(i,5)
            C_Tensor(3,1) = TDg(i,3) ; C_Tensor(3,2) = TDg(i,5) ; C_Tensor(3,3) = TDg(i,6)

            call VE_tensor(C_Tensor,Tve_Tensor)

            Tve_Tensor_tot(i,1) = Tve_Tensor(1,1)!TDg(i,1)
            Tve_Tensor_tot(i,2) = Tve_Tensor(1,2)!TDg(i,2)
            Tve_Tensor_tot(i,3) = Tve_Tensor(1,3)!TDg(i,3)
            Tve_Tensor_tot(i,4) = Tve_Tensor(2,2)!TDg(i,4)
            Tve_Tensor_tot(i,5) = Tve_Tensor(2,3)!TDg(i,5)
            Tve_Tensor_tot(i,6) = Tve_Tensor(3,3)!TDg(i,6)

            !S, Te
            S_Tensor(1,1) = TDg(i,7) ; S_Tensor(1,2) = TDg(i,8)  ; S_Tensor(1,3) = TDg(i,9)
            S_Tensor(2,1) = TDg(i,8) ; S_Tensor(2,2) = TDg(i,10) ; S_Tensor(2,3) = TDg(i,11)
            S_Tensor(3,1) = TDg(i,9) ; S_Tensor(3,2) = TDg(i,11) ; S_Tensor(3,3) = TDg(i,12)

            call Elastic_tensor(S_Tensor,Te_Tensor)

            Te_Tensor_tot(i,1) = Te_Tensor(1,1)!TDg(i,7)
            Te_Tensor_tot(i,2) = Te_Tensor(1,2)!TDg(i,8)
            Te_Tensor_tot(i,3) = Te_Tensor(1,3)!TDg(i,9)
            Te_Tensor_tot(i,4) = Te_Tensor(2,2)!TDg(i,10)
            Te_Tensor_tot(i,5) = Te_Tensor(2,3)!TDg(i,11)
            Te_Tensor_tot(i,6) = Te_Tensor(3,3)!TDg(i,12)

            !A, Ta
            A_Tensor(1,1) = TDg(i,13) ; A_Tensor(1,2) = TDg(i,14) ; A_Tensor(1,3) = TDg(i,15)
            A_Tensor(2,1) = TDg(i,14) ; A_Tensor(2,2) = TDg(i,16) ; A_Tensor(2,3) = TDg(i,17)
            A_Tensor(3,1) = TDg(i,15) ; A_Tensor(3,2) = TDg(i,17) ; A_Tensor(3,3) = TDg(i,18)

            call Orientation_tensor(A_Tensor,Ta_Tensor)

            Ta_Tensor_tot(i,1) = Ta_Tensor(1,1)!TDg(i,13)
            Ta_Tensor_tot(i,2) = Ta_Tensor(1,2)!TDg(i,14)
            Ta_Tensor_tot(i,3) = Ta_Tensor(1,3)!TDg(i,15)
            Ta_Tensor_tot(i,4) = Ta_Tensor(2,2)!TDg(i,16)
            Ta_Tensor_tot(i,5) = Ta_Tensor(2,3)!TDg(i,17)
            Ta_Tensor_tot(i,6) = Ta_Tensor(3,3)!TDg(i,18)

        end do

        write(14,'(a)')'TITLE = "FLOW"'

        variables = '"X", "Y", "Z", "P"'
        variables = trim(adjustl(variables))//', "Cxx", "Cxy", "Cxz", "Cyy", "Cyz", "Czz"'
        variables = trim(adjustl(variables))//', "Tvexx", "Tvexy", "Tvexz", "Tveyy", "Tveyz", "Tvezz"'
        variables = trim(adjustl(variables))//', "Sxx", "Sxy", "Sxz", "Syy", "Syz", "Szz"'
        variables = trim(adjustl(variables))//', "Texx", "Texy", "Texz", "Teyy", "Teyz", "Tezz"'
        variables = trim(adjustl(variables))//', "Axx", "Axy", "Axz", "Ayy", "Ayz", "Azz"'
        variables = trim(adjustl(variables))//', "Taxx", "Taxy", "Taxz", "Tayy", "Tayz", "Tazz"'
        write(14,'(2a)') 'VARIABLES=', adjustl(trim(variables))

        ! write(14,*) 'ZONE T="','TIME =',iter,'", DATAPACKING=BLOCK, N= ', NODTOL, ', E= ', NEL_3d, &
        !     ', ZONETYPE=FETETRAHEDRON, SOLUTIONTIME=', iter
        write(14,1) iter, NODTOL, NEL_3d, iter

        write(14,*)

        write(14,*) (TL(i,1),i=1,NODTOL)
        write(14,*) (TL(i,2),i=1,NODTOL)
        write(14,*) (TL(i,3),i=1,NODTOL)
        write(14,*) (TL(i,4),i=1,NODTOL)

        ! write(14,*) ( GveN*(TDg(i,1)-1.0d0), i=1,NODTOL ) !Tve,xx
        ! write(14,*) ( GveN*(TDg(i,2)-0.0d0), i=1,NODTOL ) !Tve,xy
        ! write(14,*) ( GveN*(TDg(i,3)-0.0d0), i=1,NODTOL ) !Tve,xz
        ! write(14,*) ( GveN*(TDg(i,4)-1.0d0), i=1,NODTOL ) !Tve,yy
        ! write(14,*) ( GveN*(TDg(i,5)-0.0d0), i=1,NODTOL ) !Tve,yz
        ! write(14,*) ( GveN*(TDg(i,6)-1.0d0), i=1,NODTOL ) !Tve,zz
        do ieq = 1, size(Tve_Tensor_tot,2)
            write(14,*) (TDg(i,ieq), i=1,NODTOL)
            write(14,*) (Tve_Tensor_tot(i,ieq), i=1,NODTOL)
        end do

        do ieq = 1, size(Te_Tensor_tot,2)
            write(14,*) (TDg(i,ieq+6), i=1,NODTOL)
            write(14,*) (Te_Tensor_tot(i,ieq), i=1,NODTOL)
        end do

        do ieq = 1, size(Ta_Tensor_tot,2)
            write(14,*) (TDg(i,ieq+12), i=1,NODTOL)
            write(14,*) (Ta_Tensor_tot(i,ieq), i=1,NODTOL)
        end do

        write(14,*)

        do j = 1, NEL_3d
            write(14,*)(NM_MESH(j,i),i=1,NBF_3d)
        end do

        close( 14, status='KEEP' )

        1 format('ZONE T="','TIME =', f12.4, '", STRANDID=1, DATAPACKING=BLOCK, N= ', i12, ', E= ', i12, &
            ', ZONETYPE=FETETRAHEDRON, SOLUTIONTIME=', f12.4)

    end subroutine writeTecplotFile

    !-----------------------------------------------------------------------

    subroutine writeSolution(fn, TL, TD)

        implicit none

        character(*), intent(in) :: fn
        real(8), dimension(:,:), intent(in) :: TL, TD

        integer :: inod, ieq, ierror, iel
        character(len=100) :: fn_full

        fn_full = trim(adjustl(fn))//'.DTA'
        fn_full = trim(adjustl(fn_full))

        open(12, file=fn_full, status='unknown', action='write', iostat=ierror, &
            position='rewind', form='unformatted')

        ! do ieq = 1, size(EU)
        !     write(12) EU(ieq)
        ! end do

        do inod = 1, size(TL,1)
            do ieq = 1, size(TL,2)
                write(12) TL(inod,ieq)
            end do
        end do

        do iel = 1, size(TD,1)
            do ieq = 1, size(TD,2)
                write(12) TD(iel,ieq)
            end do
        end do

        close(12)

    end subroutine writeSolution

    !-----------------------------------------------------------------------

    subroutine readSolution(fn, TL, TD)

        implicit none

        character(*), intent(in) :: fn
        real(8), dimension(:,:), intent(inout) :: TL, TD

        integer :: inod, ieq, ierr, iel
        character(len=50) :: fn_full

        fn_full = trim(adjustl(fn))//'.DTA'
        fn_full = trim(adjustl(fn_full))

        open(13, file=fn_full, status='unknown', action='read', iostat=ierr, &
            position='rewind', form='unformatted')

        ! do ieq = 1, size(EU)
        !     read(13) EU(ieq)
        ! end do

        do inod = 1, size(TL,1)
            do ieq = 1, size(TL,2)
                read(12) TL(inod,ieq)
            end do
        end do

        do iel = 1, size(TD,1)
            do ieq = 1, size(TD,2)
                read(12) TD(iel,ieq)
            end do
        end do

        close(13, status='keep')

    end subroutine readSolution

end module PostProcess_mod
