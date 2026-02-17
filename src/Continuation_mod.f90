
module Continuation_mod

    implicit none

    contains

    subroutine solveProblem()

        use Regression_mod, only: Regression
        use PHYSICAL_MODULE, only: Xparam, Xparam0, Param_Names, Param, Parametric_Factors

        implicit none

        integer :: i, j
        logical :: check

        if (Regression) then
            call initializeSolution()
            call loopRegression()
        else
            do i = 1, size(Xparam)
                do j = 1, size(Parametric_Factors)
                    check = (i > 1) .and. (j == 1)
                    ! check = (i < size(Xparam))
                    if (check) cycle
                    Param%name = Param_Names(i)
                    Param%p => Xparam(i)
                    call initializeSolution()
                    Xparam(i) = Xparam0(i)*Parametric_Factors(j)
                    write(*,*)
                    write(*,'(a,*(es12.4))') 'Parameters = ', Xparam
                    write(*,*)
                    open(15,file='Parameters.dat',position='append')
                    write(15,'(*(es14.6))') Xparam
                    close(15,status='keep')
                    call loopContinuation()
                end do
                Xparam(:) = Xparam0(:)
            end do
        end if

        write(*,'(a)') 'Problem solved! Code termination!'

    end subroutine solveProblem

    !--------------------------------------------------------------------------

    subroutine initializeSolution()

        use PHYSICAL_MODULE, only: Read_Solution
        use TIME_INTEGRATION_MODULE, only: Increment, dt, dto, dtb, dt_max, dt_min, &
            t0 => Initial_time, Total_time
        use GLOBAL_ARRAYS_MODULE, only: TL, TLo, TLb, TLp, TD, TDo, TDb, TDp, TDg
        use MESH_MODULE, only: Xm, Ym, Zm
        use PostProcess_mod, only: writeTecplotFile, readSolution
        use Tools_mod, only: getStressComponent, extrapolationLagrange

        implicit none

        integer :: i, istart, j, k
        real(8) :: L, Lo, Lb
        character(50) :: fn, str
        character(len=:), allocatable :: str1

        Increment = 0

        Total_time = t0

        dt = dt_max
        ! dt = dt_min
        dto = dt
        dtb = dto

        TLo(:,1) = Xm(:)
        TLo(:,2) = Ym(:)
        TLo(:,3) = Zm(:)
        TLo(:,4) = 0.0d0

        TDo(:,:) = 0.0d0
        istart = 1
        do i = istart, istart+5
            call getStressComponent(i, istart, j, k)
            if (j == k) then
                TDo(:,i) = 1.0d0
                TDo(:,i+6) = 1.0d0
                TDo(:,i+12) = 1.0d0/3.0d0
                TDg(:,i) = 1.0d0
                TDg(:,i+6) = 1.0d0
                TDg(:,i+12) = 1.0d0/3.0d0
            end if
        end do

        TLb(:,:) = TLo(:,:)
        TL(:,:) = TLo(:,:)
        TLp(:,:) = TLo(:,:)

        TDb(:,:) = TDo(:,:)
        TD(:,:) = TDo(:,:)
        TDp(:,:) = TDo(:,:)

        if (Read_Solution) then

            write(str,'(f12.4)') Total_time
            str1 = trim(adjustl(str))

            fn = 'SOL_'//str1
            call readSolution(fn, TL, TD)

            fn = 'SOL_'//str1//'_o'
            call readSolution(fn, TLo, TDo)

            fn = 'SOL_'//str1//'_b'
            call readSolution(fn, TLb, TDb)

            call extrapolationLagrange(t0+dt, t0-(dt+dt), t0-dt, t0, Lb, Lo, L)
            TLp(:,:) = Lb*TLb(:,:) + Lo*TLo(:,:) + L*TL(:,:)
            TDp(:,:) = Lb*TDb(:,:) + Lo*TDo(:,:) + L*TD(:,:)
            ! EUp(:) = Lb*EUb(:) + Lo*EUo(:) + L*EU(:)

            TLb(:,:) = TLo(:,:)
            TLo(:,:) = TL(:,:)
            TL(:,:) = TLp(:,:)

            TDb(:,:) = TDo(:,:)
            TDo(:,:) = TD(:,:)
            TD(:,:) = TDp(:,:)

            ! EUb(:) = EUo(:)
            ! EUo(:) = EU(:)
            ! EU(:) = EUp(:)

        end if

        call writeTecplotFile(Total_time)
        
    end subroutine initializeSolution

    !--------------------------------------------------------------------------

    subroutine loopContinuation()

        use TIME_INTEGRATION_MODULE, only: Increment, Total_time, &
            Final_time, adjust_dt
        use Tools_mod, only: writeElapsedTime, plotGraphs
        use PostProcess_mod, only: postProcess
        use Newton_mod, only: loopNewtonRapshon
        use omp_lib, only: omp_get_wtime

        implicit none

        real(8), parameter :: eps = 1.0d-8
        logical, parameter :: NR_output = .true.
        logical :: check
        real(8) :: tStart, tStart2, tEnd

        tStart = omp_get_wtime()

        timeIntegrationLoop: do

            tStart2 = omp_get_wtime()

            call printStepInfo()

            call loopNewtonRapshon(NR_output)

            call postProcess(Increment)

            if (adjust_dt) then
                call adjustTimeStep(Increment)
            end if

            call updateSolution(Increment)

            tEnd = omp_get_wtime()
            write(*,'(a)') '----------------------------------------------------------------------'
            write(*,'(a,t15,f9.4,a)') 'Step time: ', tEnd-tStart2, ' s'
            call writeElapsedTime(tStart,tEnd,'Problem time')

            if (Total_time > Final_time-eps) then
                exit timeIntegrationLoop
            end if

        end do timeIntegrationLoop

        ! check = plotGraphs(Increment)

    end subroutine loopContinuation

    !--------------------------------------------------------------------------

    subroutine loopRegression()

        use PHYSICAL_MODULE, only: Nparam
        use TIME_INTEGRATION_MODULE, only: Increment, Total_time, dt, dtb, dto, dt_max
        use Tools_mod, only: writeElapsedTime
        use Regression_mod, only: t_exp, e_exp, residualRegression, &
            jacobianRegression, solveRegressionSystem, &
            Nconstraints, residualConstraints, jacobianConstraints
        use PostProcess_mod, only: writeTecplotFile, specimenSize, writeSolution
        use Newton_mod, only: loopNewtonRapshon
        use omp_lib, only: omp_get_wtime
        use GLOBAL_ARRAYS_MODULE, only: TL, TLo, TLb, TD, TDo, TDb

        implicit none

        real(8), dimension(size(t_exp)+Nconstraints,Nparam) :: Jacobian
        real(8), dimension(size(t_exp)+Nconstraints) :: Residual
        real(8) :: tStart, tStart2, tEnd, e_eng, dt_max_loc
        integer :: s, i
        logical :: convergedRegression, data_point
        character(50) :: fn, str
        character(len=:), allocatable :: str1

        tStart = omp_get_wtime()

        s = 1   !Experimental data counter

        dt_max_loc = dt_max
        timeIntegrationLoop: do

            tStart2 = omp_get_wtime()

            call printStepInfo()

            data_point = ( abs(Total_time-t_exp(s)) < 1.0d-6 )

            call loopNewtonRapshon(data_point)

            if (data_point) then

                call writeTecplotFile(Total_time)

                e_eng = specimenSize()

                write(str,'(f12.4)') Total_time
                str1 = trim(adjustl(str))

                fn = 'out/SOL/SOL_'//str1
                call writeSolution(fn, TL, TD)

                fn = 'out/SOL/SOL_'//str1//'_o'
                call writeSolution(fn, TLo, TDo)

                fn = 'out/SOL/SOL_'//str1//'_b'
                call writeSolution(fn, TLb, TDb)

                deallocate(str1)

                Residual(s) = residualRegression(s, data_point)
                Jacobian(s,:) = jacobianRegression(s,Residual(s))

                s = s+1
                if (s > size(t_exp)) then
                    do i = 1, Nconstraints
                        Residual(s) = residualConstraints(s,data_point)
                        Jacobian(s,:) = jacobianConstraints(s,Residual(s))
                        s = s+1
                    end do
                    exit timeIntegrationLoop
                end if

            end if

            dtb = dto
            dto = dt
            if (e_exp(s) > 100.0d0) then
                dt_max_loc = 0.9999d0*dt_max_loc
            end if
            if (dt < dt_max_loc) then
                dt = 1.5d0*dt
            end if
            if (dt > dt_max_loc) dt = dt_max_loc
            dt = min(dt, t_exp(s)-Total_time)

            call updateSolution(Increment)

            tEnd = omp_get_wtime()
            write(*,'(a)') '----------------------------------------------------------------------'
            write(*,'(a,t15,f9.4,a)') 'Step time: ', tEnd-tStart2, ' s'
            call writeElapsedTime(tStart,tEnd,'Problem time')

        end do timeIntegrationLoop

        convergedRegression = solveRegressionSystem(Residual, Jacobian)
        if (.not.convergedRegression) then
            write(*,*)
            write(*,'(a)') 'Parameters updated. Re-solving the problem.'
            write(*,*)
            ! stop
            str = 'rm out/*.DAT'
            call execute_command_line(str)
            call solveProblem()
        end if

    end subroutine loopRegression

    !--------------------------------------------------------------------------

    subroutine printStepInfo()

        use TIME_INTEGRATION_MODULE, only: Increment, Total_time, dt

        implicit none

        Increment = Increment + 1

        Total_time = Total_time + dt

        write(*,'(a)')'======================================================================'
        write(*,1) Increment, Total_time, dt
        write(*,'(a)')'======================================================================'

        1 format('STEP = ', i5, 5x, ' TIME = ', f12.5, 5x, ' dt = ', f10.5)

    end subroutine printStepInfo

    !--------------------------------------------------------------------------

    subroutine adjustTimeStep(Increment)

        use TIME_INTEGRATION_MODULE, only: dt, dto, dtb, dt_min, dt_max, eps_pr
        use ELEMENTS_MODULE, only: NEQ_f
        use GLOBAL_ARRAYS_MODULE, only: TL, TLp

        implicit none

        integer, intent(in) :: Increment

        real(8) :: err_pr, trns
        real(8), dimension(NEQ_f,3) :: erind
        integer :: i

        trns = 1.0d0
        if (Increment > 3) then
            do i = 1, NEQ_f
                erind(i,1) = dot_product(TL(:,i)-TLp(:,i), TL(:,i)-TLp(:,i))
                erind(i,2) = dot_product(TL(:,i), TL(:,i))
                erind(i,3) = sqrt(erind(i,1)/erind(i,2))
            end do

            err_pr = maxval( ERIND(1:NEQ_f,3) )

            trns = sqrt(eps_pr/err_pr)

            write(*,'(a, es12.4, a, f12.6)') 'Error in prediction:', err_pr, ', trns:', trns
        end if

        dtb = dto
        dto = dt
        dt = dt*trns

        if (dt < dt_min) dt = dt_min
        if (dt > dt_max) dt = dt_max

    end subroutine adjustTimeStep

    !--------------------------------------------------------------------------

    subroutine updateSolution(Increment)

        use TIME_INTEGRATION_MODULE, only: dt, dto, dtb, Total_time
        use GLOBAL_ARRAYS_MODULE, only: TL, TLo, TLb, TLp, TD, TDo, TDb, TDp
        use Tools_mod, only: extrapolationLagrange

        implicit none

        integer, intent(in) :: Increment

        real(8) :: Lb, Lo, L

        Lb = 0.0d0 ; Lo = 0.0d0 ; L = 1.0d0
        if (Increment >= 3) then
            call extrapolationLagrange(Total_time+dt, &
                Total_time-(dtb+dto), Total_time-dto, Total_time, Lb, Lo, L)
        end if

        TLp = Lb*TLb + Lo*TLo + L*TL
        TLb = TLo
        TLo = TL
        TL  = TLp

        TDp = Lb*TDb + Lo*TDo + L*TD
        TDb = TDo
        TDo = TD
        TD  = TDp

    end subroutine updateSolution

end module Continuation_mod
