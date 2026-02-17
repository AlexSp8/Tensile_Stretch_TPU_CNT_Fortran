
module Regression_mod

    implicit none

    logical, parameter :: Regression = .true., true_stress = .false.
    !'data/Stress_neat_black.DTA', 'data/Stress_CNT_blue.DTA', 'data/Stress_CB_blue.DTA'
    !'data/Stress_ABS_hot.DTA'
    !'data/Stress_neat.DTA', 'data/Stress_CNT2.DTA'
    character(*), parameter :: data_filename = 'data/Stress_neat.DTA'

    real(8), allocatable, dimension(:) :: t_exp, e_exp, stress_exp

    integer, parameter :: Nconstraints = 0
    integer, dimension(Nconstraints), parameter :: eq_const = 0![3]

	contains

    subroutine allocateRegressionSystem()

        use Tools_mod, only: readData

        implicit none

        call readData(data_filename, t_exp, e_exp, stress_exp)

    end subroutine allocateRegressionSystem

    !---------------------------------------------------------------------

    function residualRegression(s, output) result(Res)

        use PostProcess_mod, only: stress_strain, stress_strain_true

        implicit none

        integer, intent(in) :: s
        logical, intent(in) :: output
        real(8) :: Res

        real(8) :: stress_eng
        logical :: check

        if (true_stress) then
            stress_eng = stress_strain_true(e_exp(s)/100.0d0, output)
        else
            stress_eng = stress_strain(e_exp(s)/100.0d0, output)
        end if

        Res = stress_exp(s)-stress_eng
        ! Res = log10(stress_exp(s))-log10(stress_eng)

        ! check = (data_filename == 'data/Stress_CNT_blue.DTA')
        ! check = check .and. (e_exp(s) > 110.0d0)
        ! if (check) Res = 0.0d0

        ! check = (data_filename == 'data/Stress_CNT_black.DTA')
        ! check = check .and. (s > 18)
        ! if (check) Res = 0.0d0

        ! check = (data_filename == 'data/Stress_CB_black.DTA')
        ! check = check .and. (e_exp(s) > 50.0d0)
        ! if (check) Res = 0.0d0

        ! check = (data_filename == 'data/Stress_CB_red.DTA')
        ! check = check .and. (e_exp(s) > 50.0d0)
        ! if (check) Res = 0.0d0

        ! check = (data_filename == 'data/ABS_hot.DTA')
        ! check = check .and. (e_exp(s) > 1.1d0)
        ! if (check) Res = 0.0d0

        ! check = (data_filename == 'data/ABS_cold.DTA')
        ! check = check .and. (e_exp(s) > 1.7d0)
        ! if (check) Res = 0.0d0

        ! check = (data_filename == 'data/ABS_room.DTA')
        ! check = check .and. (e_exp(s) > 1.2d0)
        ! if (check) Res = 0.0d0

        ! if (e_exp(s) > 0.99d0) Res = 1.0d4*Res
        ! if (e_exp(s) > 1.3d0) Res = 1.0d2*Res

        if (output) then
            open(14,file='eng.dat',position='append')
            write(14,'(i5,*(es18.8))') s, e_exp(s), abs(Res)
            close(14, status='keep')
        end if

    end function residualRegression

    !---------------------------------------------------------------------

    function jacobianRegression(s, Res) result(Jac)

        use PHYSICAL_MODULE, only: Xparam
        use ELEMENTS_MODULE, only: NEL_3d, NODTOL, NEQ_f, NEQ_s
        use GLOBAL_ARRAYS_MODULE, only: TL, TD
        use Newton_mod, only: loopNewtonRapshon

        implicit none

        integer, intent(in) :: s
        real(8), intent(in) :: Res
        real(8), dimension(size(Xparam)) :: Jac

        logical, parameter :: output = .false.
        real(8), parameter :: e = 1.0d-8
        real(8), dimension(size(Xparam)) :: x0
        real(8), dimension(NODTOL,NEQ_f) :: TL_initial
        real(8), dimension(NEL_3d,NEQ_s) :: TD_initial
        real(8) :: Res_temp
        integer :: j

        x0(:) = Xparam(:)

        TL_initial(:,:) = TL(:,:)
        TD_initial(:,:) = TD(:,:)

        do j = 1, size(Xparam)

            ! if (j /= getParameterID('EpN')) then
                Xparam(j) = Xparam(j) + e
            ! end if

            call loopNewtonRapshon(output)

            Res_temp = residualRegression(s, output)

            Jac(j) = (Res_temp - Res)/e

            Xparam(j) = x0(j)

            TL(:,:) = TL_initial(:,:)
            TD(:,:) = TD_initial(:,:)

        end do

    end function jacobianRegression

    !------------------------------------------------------------------------------

    logical function solveRegressionSystem(Residual, Jacobian) result(convergedRegression)

        use PHYSICAL_MODULE, only: Xparam
        use Tools_mod, only: getParameterID, plotGraphs

        implicit none

        real(8), dimension(:), intent(in) :: Residual
        real(8), dimension(:,:), intent(in) :: Jacobian

        integer, parameter :: N = size(Xparam)
        real(8), parameter :: NR_Tolerance = 1.0d-6
        integer, save :: iSolves = 0
        real(8), dimension(N) :: x0, JTR, sol
        real(8), dimension(N,N) :: JTJ
        real(8), dimension(N,size(Residual)) :: JT
        integer, dimension(N) :: ipiv

        integer :: i, info
        real(8) :: relax, JTR_norm, Cor_norm, R_norm
        logical :: check

        convergedRegression = .false.

        relax = 1.0d0

        JT = transpose(Jacobian)
        JTJ = matmul(JT,Jacobian)
        JTR = matmul(JT,Residual)

        JTR_norm = norm2(JTR)
        R_norm = norm2(Residual)

        check = (JTR_norm < NR_Tolerance) .or. (R_norm < NR_Tolerance)
        if (check) then
            write(*,'(a)')'Non linear regression converged!'
            write(*,'(a,es12.4)')'JTR norm:        ', JTR_norm
            write(*,'(a,es12.4)')'R norm:          ', R_norm
            convergedRegression = .true.
        end if

        do i = 1, N
            JTJ(i,i) = JTJ(i,i) + JTR_norm
        end do

        !perform LU decomposition  
        call DGETRF(N, N, JTJ, N, ipiv, info)

        !perform back-substitution
        call DGETRS('N', N, 1, JTJ, N, ipiv, JTR, N, info)

        !Copy "b" to "sol"
        sol(:) = -JTR(:)

        Cor_norm = norm2(sol)

        x0(:) = Xparam(:)

        Xparam(:) = Xparam(:) + relax*sol(:)

        ! i = getParameterID('TdN')
        ! Xparam(i) = x0(i)
        
        ! i = getParameterID('GveN')
        ! Xparam(i) = x0(i)

        ! i = getParameterID('GeN')
        ! Xparam(i) = x0(i)

        ! i = getParameterID('JmN')
        ! Xparam(i) = x0(i)

        i = getParameterID('EpN')
        Xparam(i) = x0(i)

        ! i = getParameterID('TaN')
        ! Xparam(i) = x0(i)

        if (mod(iSolves,1) == 0) then
            check = plotGraphs(iSolves)
        end if
        iSolves = iSolves+1

        write(*,'(a)') '----------------------------------------------------------------------'
        write(*,'(a,es12.4)')'JTR norm        ', JTR_norm
        write(*,'(a,es12.4)')'R norm          ', R_norm
        write(*,'(a,es12.4)')'Correction norm ', Cor_norm
        write(*,'(a,es12.4)')'Max Res         ', maxval(abs(Residual))

        open(15,file='Parameters.dat',position='append')
        write(15,'(*(es14.6))') Xparam, JTR_norm, R_norm
        close(15,status='keep')

        write(*,*)
        write(*,'(a,*(es14.6))') 'x = ', Xparam, JTR_norm, R_norm

    end function solveRegressionSystem

    !------------------------------------------------------------------------------

    function residualConstraints(s, output) result(Res)
        
        use PHYSICAL_MODULE, only: x0 => Xparam0, x => Xparam

        implicit none

        integer, intent(in) :: s
        logical, intent(in) :: output
        real(8) :: Res

        integer :: iC, Nexp, jeq_C
        real(8) :: penalty, f, ub, lb

        Nexp = size(t_exp)
        iC = s - Nexp

        jeq_C = eq_const(iC)

        !boundary constraints: up to 20% difference
        ub = 2.0d0!1.2d0
        lb = 1.0d0!0.8d0
        penalty = 1.0d1

        !Lower bound
        ! f = 1.0d0/(lb*x0(jeq_C)-x(jeq_C))
        f = max( 0.0d0, lb*x0(jeq_C)-x(jeq_C) )
        Res = penalty*(f**2)

        if (output .and. f > 0.0d0) then
            write(*,'(i3, 2es12.4)') jeq_C, -f, Res
        end if

        ! !Upper bound
        ! ! f = 1.0d0/(x(jeq_C)-ub*x0(jeq_C))
        ! f = max( 0.0d0, x(jeq_C)-ub*x0(jeq_C) )
        ! Res = penalty*(f**2)

        ! if (output .and. f > 0.0d0) then
        !     write(*,'(2i3, f12.4, es12.4)') jeq_C, f, Res
        ! end if

    end function residualConstraints

    !---------------------------------------------------------------------

    function jacobianConstraints(s, Res) result(Jac)

        use PHYSICAL_MODULE, only: x => Xparam

        implicit none

        integer, intent(in) :: s
        real(8), intent(in) :: Res
        real(8), dimension(size(x)) :: Jac

        logical, parameter :: output = .false.
        real(8), parameter :: e = 1.0d-8
        real(8), dimension(size(x)) :: x0
        real(8) :: Res_temp
        integer :: j

        x0(:) = x(:)

        do j = 1, size(x)

            ! if (j /= getParameterID('EpN')) then
                x(j) = x(j) + e
            ! end if

            Res_temp = residualConstraints(s, output)

            Jac(j) = (Res_temp - Res)/e

            x(j) = x0(j)

        end do

    end function jacobianConstraints

end module Regression_mod
