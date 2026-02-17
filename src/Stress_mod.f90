
module Stress_mod

    implicit none

    contains
    
    subroutine DOMI_STRESS(iel, TEMP_TL, TEMP_TD)

        use TIME_INTEGRATION_MODULE, only: Increment, dt
        use ELEMENTS_MODULE, only: NBF_3d, NEQ_f, NEQ_s
        use GAUSS_MODULE, only: LINEAR_TETRAHEDRON
        use ENUMERATION_MODULE, only: NM_MESH
        use MESH_MODULE, only: Xm, Ym, Zm
        use GLOBAL_ARRAYS_MODULE, only: TLo, TLb
        use Tools_mod, only: basis3D

        implicit none

        integer, intent(in) :: iel
        real(8), dimension(NBF_3d, NEQ_f), intent(in)  :: TEMP_TL
        real(8), dimension(NEQ_s), intent(out) :: TEMP_TD

        integer :: inod, gnod, i
        real(8) :: CJAC
        real(8), dimension(3) :: Xi_gs, dt_term
        real(8), dimension(NBF_3d) :: bfn
        real(8), dimension(3,3) :: GV, dXidci
        real(8), dimension(NBF_3d,3) :: Xi_el, dfdci, dfdxi

        !BASIS FUNCTIONS AT THE CENTROID OF THE ELEMENT
        call LINEAR_TETRAHEDRON( bfn, dfdci(:,1), dfdci(:,2), dfdci(:,3), &
            1.0d0/4.0d0, 1.0d0/4.0d0, 1.0d0/4.0d0 )

        Xi_el(:,:) = TEMP_TL(:,1:3)
        ! Xi_el(:,1) = Xm(NM_MESH(iel,:))
        ! Xi_el(:,2) = Ym(NM_MESH(iel,:))
        ! Xi_el(:,3) = Zm(NM_MESH(iel,:))

        call basis3D(Xi_el, bfn, dfdci, dfdxi, Xi_gs, dXidci, CJAC)

        GV(:,:) = 0.0d0
        do inod = 1, NBF_3d

            gnod = NM_MESH(iel,inod)

            if (Increment == 1) then
                dt_term(:) = (TEMP_TL(inod,1:3)-TLo(gnod,1:3))/dt
            else
                dt_term(:) = (3.0d0*TEMP_TL(inod,1:3)-4.0d0*TLo(gnod,1:3)+TLb(gnod,1:3))/(2.0d0*dt)
            end if
            do i = 1, 3
                GV(:,i) = GV(:,i) + dt_term(i)*dfdxi(inod,:)
            end do

        end do

        call ELEMENTAL_STRESS( iel, TEMP_TD, GV )

    end subroutine DOMI_STRESS

    !---------------------------------------------------------------------

    subroutine ELEMENTAL_STRESS( iel, X_stress, GV )

        use ELEMENTS_MODULE, only: NEQ_s
        use NRAPSHON_MODULE, only: ERROR_STR, MNR_eps
        use GLOBAL_ARRAYS_MODULE, only: TDo, TDb
        use BOUNDARY_ENUMERATION_MODULE, only: NBE, BND_Symmetry, Symmetry

        implicit none

        integer, intent(in) :: iel
        real(8), dimension(:), intent(inout) :: X_stress
        real(8), dimension(:,:), intent(in) :: GV

        integer, parameter :: N = NEQ_s
        integer, parameter :: maxiter = 10000
        integer, dimension(N) :: ipvt
        integer, dimension(6) :: ieqs
        integer :: iter, info, ibnd, imax, i, ieq, jeq
        real(8), dimension(N,N) :: Jacobian
        real(8), dimension(N) :: Residual, WORK
        real(8), dimension(N) :: Xo_stress, Xb_stress, X0_stress
        real(8) :: Res_norm, Cor_norm, COND, xf
        character(3) :: Flag_NR

        X0_stress(:) = X_stress(:)
        Xo_stress = TDo(iel,:)
        Xb_stress = TDb(iel,:)

        iter = 0
        Flag_NR = 'NRP'
        xf = 1.0d0
        imax = maxiter
        LOOP_NEWTONRAPSHON: do

            iter = iter + 1

            if (iter > imax) then
                write(*,'(a,t28,2i10,3es12.4)') 'Stress NR Max Iter at iel:', iel, imax, &
                    xf, Res_norm, Cor_norm
                ! stop
                X_stress(:) = X0_stress(:)
                xf = xf/2.0d0
                iter = 0
                Flag_NR = 'NRP'
                imax = 2*imax
                cycle
            end if

            if (xf < 1.0d-2) exit LOOP_NEWTONRAPSHON

            call STRESS_Residual( X_stress, Xo_stress, Xb_stress, GV, Residual )

            Res_norm = dot_product(Residual, Residual)
            Res_norm = sqrt(Res_norm)

            if ( Res_norm<ERROR_STR ) then
                    ! if (iel==1) print*, iter, 'Res_norm = ', Res_norm
                exit LOOP_NEWTONRAPSHON
            end if

            if (Flag_NR == 'NRP') then
                call STRESS_Jacobian( X_stress, Xo_stress, Xb_stress, GV, Residual, Jacobian )
            end if

            !Y = 0: Cxy, Cyz, Sxy, Syz, Axy, Ayz
            do i = 1, size(BND_Symmetry)
                ibnd = BND_Symmetry(i)

                select case (Symmetry(i))

                case ('X')
                    !xy(2,8,14),xz(3,9,15)
                    ieqs = [2, 3, 8, 9, 14, 15]
                case ('Y')
                    !xy(2,8,14),yz(5,11,17)
                    ieqs = [2, 5, 8, 11, 14, 17]
                case ('Z')
                    !xz(3,9,15),yz(5,11,17)
                    ieqs = [3, 5, 9, 11, 15, 17]
                case default
                    write(*,'(a)') 'Wrong Symmetry condition in ELEMENTAL_STRESS!'
                    stop

                end select

                if (NBE(iel)%KOE(ibnd)) then
                    do ieq = 1, size(ieqs)
                        jeq = ieqs(ieq)
                        Residual(jeq) = X_stress(jeq) - 0.0d0
                        Jacobian(jeq,:) = 0.0d0 ; Jacobian(jeq,jeq) = 1.0d0
                    end do
                    ! !Cxy
                    ! Residual(2) = X_stress(2) - 0.0d0
                    ! Jacobian(2,:) = 0.0d0 ; Jacobian(2,2) = 1.0d0
                    ! !Cyz
                    ! Residual(5) = X_stress(5) - 0.0d0
                    ! Jacobian(5,:) = 0.0d0 ; Jacobian(5,5) = 1.0d0
                    ! !Sxy
                    ! Residual(8) = X_stress(8) - 0.0d0
                    ! Jacobian(8,:) = 0.0d0 ; Jacobian(8,8) = 1.0d0
                    ! !Syz
                    ! Residual(11) = X_stress(11) - 0.0d0
                    ! Jacobian(11,:) = 0.0d0 ; Jacobian(11,11) = 1.0d0

                end if
            end do

            ! call DECOMP( N, N, Jacobian, COND, ipvt, WORK )
            ! call SOLVE( N, N, Jacobian, Residual, ipvt )

            !perform LU decomposition  
            call DGETRF(N, N, Jacobian, N, ipvt, info)
            !perform back-substitution
            call DGETRS('N', N, 1, Jacobian, N, ipvt, Residual, N, info)

            X_stress(:) = X_stress(:) - xf*Residual(:)

            Cor_norm = dot_product(Residual,Residual)
            Cor_norm = sqrt(Cor_norm)

            if ( Cor_norm<ERROR_STR ) then
                    ! if (iel==1) print*, iter, 'Cor_Norm = ', Cor_norm
                exit LOOP_NEWTONRAPSHON
            end if

            if (Cor_norm < MNR_eps) then
                Flag_NR = 'MNR'
            end if

        end do LOOP_NEWTONRAPSHON

    end subroutine ELEMENTAL_STRESS

    !---------------------------------------------------------------------

    subroutine STRESS_Residual( X, Xo, Xb, GV, Residual )

        use TIME_INTEGRATION_MODULE, only: Increment, dt
        use PHYSICAL_MODULE, only: VE_model, a_UCD, Xparam, ThetaN
        use ConstitutiveModels_mod
        use Tools_mod, only: traceTensor, getParameterID, doubleDotProduct!, inverseTensor3D

        implicit none

        real(8), dimension(:), intent(in) :: X, Xo, Xb
        real(8), dimension(:,:), intent(in) :: GV
        real(8), dimension(:), intent(out) :: Residual

        integer, parameter, dimension(3,3) :: I_Tensor = [1,0,0,0,1,0,0,0,1]
        real(8), dimension(3,3) :: dCdt, UCD, C_Tensor, Co_Tensor, Cb_Tensor, GVT
        real(8), dimension(3,3) :: C2_Tensor, W_Tensor, D_Tensor, G_dot
        real(8), dimension(3) :: f
        real(8), dimension(3,3) :: dSdt, UCD_S, S_Tensor, So_Tensor, Sb_Tensor, S_term1!, S_term2, SI_Tensor, S_Conv
        real(8), dimension(3,3) :: dAdt, UCD_A, A_Tensor, Ao_Tensor, Ab_Tensor, A_term1, A_term2, A_term3
        real(8) :: divV, TdN, TaN
        integer :: i, j, k, istart

        !C tensor evolution
        TdN = Xparam(getParameterID('TdN'))
        
        C_Tensor(1,1) = X(1) ; C_Tensor(1,2) = X(2) ; C_Tensor(1,3) = X(3)
        C_Tensor(2,1) = X(2) ; C_Tensor(2,2) = X(4) ; C_Tensor(2,3) = X(5)
        C_Tensor(3,1) = X(3) ; C_Tensor(3,2) = X(5) ; C_Tensor(3,3) = X(6)

        Co_Tensor(1,1) = Xo(1) ; Co_Tensor(1,2) = Xo(2) ; Co_Tensor(1,3) = Xo(3)
        Co_Tensor(2,1) = Xo(2) ; Co_Tensor(2,2) = Xo(4) ; Co_Tensor(2,3) = Xo(5)
        Co_Tensor(3,1) = Xo(3) ; Co_Tensor(3,2) = Xo(5) ; Co_Tensor(3,3) = Xo(6)

        Cb_Tensor(1,1) = Xb(1) ; Cb_Tensor(1,2) = Xb(2) ; Cb_Tensor(1,3) = Xb(3)
        Cb_Tensor(2,1) = Xb(2) ; Cb_Tensor(2,2) = Xb(4) ; Cb_Tensor(2,3) = Xb(5)
        Cb_Tensor(3,1) = Xb(3) ; Cb_Tensor(3,2) = Xb(5) ; Cb_Tensor(3,3) = Xb(6)

        if (Increment == 1) then
            dCdt(:,:) = (C_Tensor(:,:)-Co_Tensor(:,:))/dt
        else
            dCdt(:,:) = (3.0d0*C_Tensor(:,:)-4.0d0*Co_Tensor(:,:)+Cb_Tensor(:,:))/(2.0d0*dt)
        end if

        GVT(:,:) = transpose(GV)

        ! UCD(:,:) = dCdt(:,:) - matmul(GVT,C_Tensor) - matmul(C_Tensor,GV)

        W_Tensor(:,:) = 0.5d0*(GVT(:,:)-GV(:,:))
        D_Tensor(:,:) = 0.5d0*(GV(:,:)+GVT(:,:))
        UCD(:,:) = dCdt(:,:) - matmul(W_Tensor,C_Tensor) + matmul(C_Tensor,W_Tensor) &
                    - a_UCD*(matmul(D_Tensor,C_Tensor)+matmul(C_Tensor,D_Tensor))

        select case (VE_model)
        case ('LPTT')
            call VE_model_LPTT(C_Tensor,f)
        case ('ePTT')
            call VE_model_ePTT(C_Tensor,f)
        case ('Giesekus')
            call VE_model_Giesekus(f)
        case ('FENE-CR')
            call VE_model_FENE_CR(C_Tensor,f)
        case ('FENE-P')
            call VE_model_FENE_P(C_Tensor,f)
        case default
            write(*,'(a)') 'Wrong VE_model in STRESS_Residual!'
            stop
        end select

        C2_Tensor(:,:) = matmul(C_Tensor,C_Tensor)
        ! istart = 1
        ! do i = 1, size(Residual)
        !     call getStressComponent(i, istart, j, k)
        !     Residual(i) = TdN*UCD(j,k) + f(3)*C2_Tensor(j,k) + f(2)*C_Tensor(j,k) - f(1)*I_Tensor(j,k)
        ! end do
        Residual(1) = TdN*UCD(1,1) + f(3)*C2_Tensor(1,1) + f(2)*C_Tensor(1,1) - f(1)*I_Tensor(1,1)
        Residual(2) = TdN*UCD(1,2) + f(3)*C2_Tensor(1,2) + f(2)*C_Tensor(1,2) - f(1)*I_Tensor(1,2)
        Residual(3) = TdN*UCD(1,3) + f(3)*C2_Tensor(1,3) + f(2)*C_Tensor(1,3) - f(1)*I_Tensor(1,3)
        Residual(4) = TdN*UCD(2,2) + f(3)*C2_Tensor(2,2) + f(2)*C_Tensor(2,2) - f(1)*I_Tensor(2,2)
        Residual(5) = TdN*UCD(2,3) + f(3)*C2_Tensor(2,3) + f(2)*C_Tensor(2,3) - f(1)*I_Tensor(2,3)
        Residual(6) = TdN*UCD(3,3) + f(3)*C2_Tensor(3,3) + f(2)*C_Tensor(3,3) - f(1)*I_Tensor(3,3)

        !S tensor evolution
        S_Tensor(1,1) = X(7) ; S_Tensor(1,2) = X(8)  ; S_Tensor(1,3) = X(9)
        S_Tensor(2,1) = X(8) ; S_Tensor(2,2) = X(10) ; S_Tensor(2,3) = X(11)
        S_Tensor(3,1) = X(9) ; S_Tensor(3,2) = X(11) ; S_Tensor(3,3) = X(12)

        So_Tensor(1,1) = Xo(7) ; So_Tensor(1,2) = Xo(8)  ; So_Tensor(1,3) = Xo(9)
        So_Tensor(2,1) = Xo(8) ; So_Tensor(2,2) = Xo(10) ; So_Tensor(2,3) = Xo(11)
        So_Tensor(3,1) = Xo(9) ; So_Tensor(3,2) = Xo(11) ; So_Tensor(3,3) = Xo(12)

        Sb_Tensor(1,1) = Xb(7) ; Sb_Tensor(1,2) = Xb(8)  ; Sb_Tensor(1,3) = Xb(9)
        Sb_Tensor(2,1) = Xb(8) ; Sb_Tensor(2,2) = Xb(10) ; Sb_Tensor(2,3) = Xb(11)
        Sb_Tensor(3,1) = Xb(9) ; Sb_Tensor(3,2) = Xb(11) ; Sb_Tensor(3,3) = Xb(12)

        if (Increment == 1) then
            dSdt(:,:) = (S_Tensor(:,:)-So_Tensor(:,:))/dt
        else
            dSdt(:,:) = (3.0d0*S_Tensor(:,:)-4.0d0*So_Tensor(:,:)+Sb_Tensor(:,:))/(2.0d0*dt)
        end if

        ! UCD_S(:,:) = dSdt(:,:) - matmul(GVT,S_Tensor) - matmul(S_Tensor,GV)
        UCD_S(:,:) = dSdt(:,:) - matmul(W_Tensor,S_Tensor) + matmul(S_Tensor,W_Tensor) &
                    - a_UCD*(matmul(D_Tensor,S_Tensor)+matmul(S_Tensor,D_Tensor))

        divV = traceTensor(GV)
        S_term1(:,:) = (2.0d0/3.0d0)*divV*S_Tensor(:,:)

        ! call inverseTensor3D(S_Tensor, SI_Tensor)
        ! S_Conv(:,:) = 0.0d0
        ! do k = 1, 3
        !     S_Conv(:,:) = S_Conv(:,:) + V_f(k)*dSdXi(:,:,k)
        ! end do
        ! S_term2(:,:) = doubleDotProduct(SI,S_Conv)*S_Tensor(:,:)

        Residual(7)  = UCD_S(1,1) + S_term1(1,1) !+ S_term2(1,1)
        Residual(8)  = UCD_S(1,2) + S_term1(1,2) !+ S_term2(1,2)
        Residual(9)  = UCD_S(1,3) + S_term1(1,3) !+ S_term2(1,3)
        Residual(10) = UCD_S(2,2) + S_term1(2,2) !+ S_term2(2,2)
        Residual(11) = UCD_S(2,3) + S_term1(2,3) !+ S_term2(2,3)
        Residual(12) = UCD_S(3,3) + S_term1(3,3) !+ S_term2(3,3)

        !A tensor evolution
        TaN = Xparam(getParameterID('TaN'))

        A_Tensor(1,1) = X(13) ; A_Tensor(1,2) = X(14) ; A_Tensor(1,3) = X(15)
        A_Tensor(2,1) = X(14) ; A_Tensor(2,2) = X(16) ; A_Tensor(2,3) = X(17)
        A_Tensor(3,1) = X(15) ; A_Tensor(3,2) = X(17) ; A_Tensor(3,3) = X(18)

        Ao_Tensor(1,1) = Xo(13) ; Ao_Tensor(1,2) = Xo(14) ; Ao_Tensor(1,3) = Xo(15)
        Ao_Tensor(2,1) = Xo(14) ; Ao_Tensor(2,2) = Xo(16) ; Ao_Tensor(2,3) = Xo(17)
        Ao_Tensor(3,1) = Xo(15) ; Ao_Tensor(3,2) = Xo(17) ; Ao_Tensor(3,3) = Xo(18)

        Ab_Tensor(1,1) = Xb(13) ; Ab_Tensor(1,2) = Xb(14) ; Ab_Tensor(1,3) = Xb(15)
        Ab_Tensor(2,1) = Xb(14) ; Ab_Tensor(2,2) = Xb(16) ; Ab_Tensor(2,3) = Xb(17)
        Ab_Tensor(3,1) = Xb(15) ; Ab_Tensor(3,2) = Xb(17) ; Ab_Tensor(3,3) = Xb(18)

        if (Increment == 1) then
            dAdt(:,:) = (A_Tensor(:,:)-Ao_Tensor(:,:))/dt
        else
            dAdt(:,:) = (3.0d0*A_Tensor(:,:)-4.0d0*Ao_Tensor(:,:)+Ab_Tensor(:,:))/(2.0d0*dt)
        end if

        ! UCD_A(:,:) = dAdt(:,:) - matmul(GVT,A_Tensor) - matmul(A_Tensor,GV)
        UCD_A(:,:) = dAdt(:,:) - matmul(W_Tensor,A_Tensor) + matmul(A_Tensor,W_Tensor) &
                    - a_UCD*(matmul(D_Tensor,A_Tensor)+matmul(A_Tensor,D_Tensor))

        A_term1(:,:) = doubleDotProduct(A_Tensor,GV)*A_Tensor(:,:)
        ! A_term1(:,:) = matmul(A_Tensor(:,:),matmul(A_Tensor,GV))
        A_term1(:,:) = 2.0d0*ThetaN*A_term1(:,:)

        G_dot(:,:) = GV(:,:) + GVT(:,:)
        A_term2(:,:) = matmul(A_Tensor,G_dot)+matmul(G_dot,A_Tensor)
        A_term2(:,:) = (ThetaN-1.0d0)*A_term2(:,:)/2.0d0

        A_term3(:,:) = A_Tensor(:,:) - I_Tensor(:,:)/3.0d0
        A_term3(:,:) = A_term3(:,:)/TaN

        Residual(13) = UCD_A(1,1) + A_term1(1,1) + A_term2(1,1) + A_term3(1,1)
        Residual(14) = UCD_A(1,2) + A_term1(1,2) + A_term2(1,2) + A_term3(1,2)
        Residual(15) = UCD_A(1,3) + A_term1(1,3) + A_term2(1,3) + A_term3(1,3)
        Residual(16) = UCD_A(2,2) + A_term1(2,2) + A_term2(2,2) + A_term3(2,2)
        Residual(17) = UCD_A(2,3) + A_term1(2,3) + A_term2(2,3) + A_term3(2,3)
        Residual(18) = UCD_A(3,3) + A_term1(3,3) + A_term2(3,3) + A_term3(3,3)

    end subroutine STRESS_Residual

    !-----------------------------------------------------------------------

    subroutine STRESS_Jacobian( X, Xo, Xb, GradU, Residual, Jacobian )

        use Tools_mod, only: perturbVariable

        implicit none

        real(8), dimension(:), intent(inout) :: X
        real(8), dimension(:), intent(in) :: Xo, Xb
        real(8), dimension(:,:), intent(in) :: GradU
        real(8), dimension(:), intent(in) :: Residual
        real(8), dimension(:,:), intent(out) :: Jacobian

        integer :: i
        real(8) :: EPS_JAC
        real(8), dimension(size(Residual)) :: Resf, Resb

        ! do i = 1, size(X)

        !     EPS_JAC = perturbVariable( X(i) )

        !     X(i) = X(i) + EPS_JAC

        !     call STRESS_Residual( X, Xo, Xb, GradU, Resf )

        !     X(i) = X(i) - 2.0d0*EPS_JAC

        !     call STRESS_Residual( X, Xo, Xb, GradU, Resb )

        !     Jacobian(:,i) = (Resf-Resb)/(2.0d0*EPS_JAC)

        !     X(i) = X(i) + EPS_JAC

        ! end do

        do i = 1, size(X)

            EPS_JAC = perturbVariable( X(i) )

            X(i) = X(i) + EPS_JAC

            call STRESS_Residual( X, Xo, Xb, GradU, Resf )

            Jacobian(:,i) = (Resf-Residual)/EPS_JAC

            X(i) = X(i) - EPS_JAC

        end do

    end subroutine STRESS_Jacobian

    !-----------------------------------------------------------------------

    subroutine VE_tensor(C_Tensor, Tve_Tensor)

        use PHYSICAL_MODULE, only: VE_model, Xparam
        use ConstitutiveModels_mod, only: VE_tensor_regular, &
            VE_tensor_FENE_CR, VE_tensor_FENE_P
        use Tools_mod, only: getParameterID

        implicit none

        real(8), dimension(:,:), intent(in) :: C_Tensor
        real(8), dimension(:,:), intent(out) :: Tve_Tensor

        integer, parameter, dimension(3,3) :: I_Tensor = [1,0,0,0,1,0,0,0,1]
        real(8), dimension(2) :: g
        real(8) :: GveN

        select case (VE_model)
        case ('LPTT', 'ePTT', 'Giesekus')
            call VE_tensor_regular(g)
        case ('FENE-CR')
            call VE_tensor_FENE_CR(C_Tensor,g)
        case ('FENE-P')
            call VE_tensor_FENE_P(C_Tensor,g)
        case default
            write(*,'(a)') 'Wrong VE_model in VE_tensor'
            stop
        end select

        GveN = Xparam(getParameterID('GveN'))

        Tve_Tensor(:,:) = GveN*(g(2)*C_Tensor(:,:) - g(1)*I_Tensor(:,:))

    end subroutine VE_tensor

    !-----------------------------------------------------------------------

    subroutine Elastic_tensor(S_Tensor, Te_Tensor)

        use PHYSICAL_MODULE, only: Elastic_model, Ndim, Xparam
        use Tools_mod, only: traceTensor, determinant2D, getParameterID
        use ConstitutiveModels_mod, only: Elastic_tensor_Neo_Hookean, Elastic_tensor_Mooney, &
            Elastic_tensor_Gent_I, Elastic_tensor_Gent_Thomas, Elastic_tensor_Gent_II

        implicit none

        real(8), dimension(:,:), intent(in) :: S_Tensor
        real(8), dimension(:,:), intent(out) :: Te_Tensor

        integer, parameter, dimension(3,3) :: I_Tensor = [1,0,0,0,1,0,0,0,1]
        real(8), dimension(3,3) :: S2_Tensor
        real(8), dimension(3) :: dwi, Ii
        real(8) :: trace, trace2, W, GeN, JmN, BetaN, BetaN2

        GeN = Xparam(getParameterID('GeN'))

        trace = traceTensor(S_Tensor)
        S2_Tensor(:,:) = matmul(S_Tensor,S_Tensor)
        trace2 = traceTensor(S2_Tensor)

        Ii(1) = trace
        Ii(2) = 0.5d0*(trace**2-trace2)
        Ii(3) = determinant2D(S_Tensor)

        select case (Elastic_model)
        case ('Neo-Hookean')
            W = Ii(1)-3.0d0
            call Elastic_tensor_Neo_Hookean(Ii, dwi)
        case ('Mooney')
            BetaN = Xparam(getParameterID('BetaN'))
            W = Ii(1)-3.0d0 + BetaN*(Ii(2)-3.0d0)
            call Elastic_tensor_Mooney(Ii, dwi)
        case ('Gent-Thomas')
            BetaN = Xparam(getParameterID('BetaN'))
            W = Ii(1)-3.0d0 + BetaN*log(Ii(2)/3.0d0)
            call Elastic_tensor_Gent_Thomas(Ii, dwi)
        case ('Gent-I')
            JmN = Xparam(getParameterID('JmN'))
            W = -JmN*log(1.0d0-(Ii(1)-3.0d0)/JmN)
            call Elastic_tensor_Gent_I(Ii, dwi)
        case ('Gent-II')
            JmN = Xparam(getParameterID('JmN'))
            BetaN2 = Xparam(getParameterID('BetaN2'))
            W = -JmN*log(1.0d0-(Ii(1)-3.0d0)/JmN) + BetaN2*log(Ii(2)/3.0d0)
            call Elastic_tensor_Gent_II(Ii, dwi)
        case default
            write(*,'(a)') 'Wrong Elastic_model in Elastic_tensor'
            stop
        end select

        Te_Tensor(:,:) = GeN*( (dwi(1)+Ii(1)*dwi(2))*S_Tensor(:,:) - dwi(2)*S2_Tensor(:,:) &
                                - (1.0d0/Ndim)*(Ii(1)*dwi(1)+2.0d0*Ii(2)*dwi(2))*I_Tensor(:,:) )

    end subroutine Elastic_tensor

    !-----------------------------------------------------------------------

    subroutine Orientation_tensor(A_Tensor, Ta_Tensor)

        use PHYSICAL_MODULE, only: TempN, KbN, PhiN, LpN, DpN, PI, ThetaN

        implicit none

        real(8), dimension(:,:), intent(in) :: A_Tensor
        real(8), dimension(:,:), intent(out) :: Ta_Tensor
        integer, parameter, dimension(3,3) :: I_Tensor = [1,0,0,0,1,0,0,0,1]

        real(8), parameter :: na0 = 6.0d0/(PI*LpN*DpN**2) !in mm^-3
        real(8) :: factor, na
        real(8), dimension(3) :: temp
        integer :: i

        na = na0*PhiN*1.0d9 !in m^-3

        factor = 3.0d0*na*KbN*TempN*ThetaN !in kg m^-1 s^-2 -> Pa
        factor = 1.0d-6*factor !in MPa

        Ta_Tensor(:,:) = factor*(A_Tensor(:,:)-I_Tensor(:,:)/3.0d0)

        ! do i = 1, 3
        !     temp = A_Tensor(i,:)
        !     write(*,'(*(es24.12))') temp
        ! end do
        ! do i = 1, 3
        !     temp = Ta_Tensor(i,:)
        !     write(*,'(*(es24.12))') temp
        ! end do
        ! pause

    end subroutine Orientation_tensor

end module Stress_mod
