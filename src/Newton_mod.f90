
module Newton_mod

    implicit none

    contains

    subroutine loopNewtonRapshon(output)

        use TIME_INTEGRATION_MODULE, only: Increment, dt, dto
        use ELEMENTS_MODULE, only: NEL_3d, NODTOL, NEQ_f, NUNKNOWNS_f
        use NRAPSHON_MODULE, only: NITER, ERROR_NR, MNR_eps, iter0, Cor_norm0
        use FLOW_ARRAYS_MODULE, only: B_f, S_f
        use GLOBAL_ARRAYS_MODULE, only: TL, TLo, TLp
        use CSR_STORAGE, only: PT_f, MAXFCT_f, MNUM_f, MTYPE_f, PHASE_f, N_f, A_f, IA_f, CA_f,&
            IDUM_f, NRHS_f, IPARM_f, MSGLVL_f, DDUM_f, DDUM_f, ERROR_f
        use OMP_PARALLEL, only: NTHREADS
        use omp_lib, only: omp_get_wtime
        use Equations_mod, only: FLOW_EQUATIONS, DIRICHLET_BOUNDARY_CONDITIONS, STRESS_STORAGE

        implicit none

        logical, intent(in) :: output

        integer :: i, j, iel, k, xNITER, ITER_f
        real(8) :: RES_NORM, COR_NORM_OLD_f, COR_NORM_NEW_f, xERROR_NR, xF
        logical :: EMERGENCY, LMSG, check
        real(8) :: tStart, tEnd, tStart2, tStart3
        character(len=3) :: FLAG_NR

        tStart = omp_get_wtime()

        EMERGENCY = .false.
        xF = 1.0d0

        300 continue

        if (EMERGENCY) then
            xNITER = 50*NITER
            xF = 0.5d0*xF
            if (xF<1.0d-2) then
                if (output) write(*,'(a)') 'VERY SMALL RELAXATION FACTOR, GENERAL STOP!'
                stop
            end if
            xERROR_NR = 1.0d0*ERROR_NR
            TL(:,:) = TLp(:,:)
            if (output) write(*,*) 'xF=', xF
        else
            xNITER = NITER
            xF = 1.0d0
            xERROR_NR = ERROR_NR
        end if

        ! if (Increment == 1) then
        !     xF = 0.5d0
        !     xNITER = 100*NITER
        ! end if

        if (Increment == 1) then
            iter0 = 0
            Cor_Norm0 = 1.0d0
        end if

        ITER_f = 0

        FLAG_NR='NRP'
        check = (Cor_Norm0 < MNR_eps/10.0d0)
        check = check .and. (iter0 < 10)
        check = check .and. (abs(dt-dto) < 1.0d-6)
        ! check = check .and. (Increment > 10)
        if (check) FLAG_NR='MNR'

        RES_NORM = 1.0d0
        COR_NORM_OLD_f = 1.0d0
        COR_NORM_NEW_f = 1.0d0

        LOOP_loopNewtonRapshon: do

            check = (COR_NORM_NEW_f < xERROR_NR)
            ! check = check .and. (RES_NORM < xERROR_NR)
            if (check) then
                iter0 = ITER_f
                exit
            end if

            tStart2 = omp_get_wtime()

            ITER_f = ITER_f + 1

            ! CHECK_EMERGENCY: if (EMERGENCY) then
            !     FLAG_NR='NRP'
            ! end if CHECK_EMERGENCY

            B_f = 0.0d0
            S_f = 0.0d0
            if (FLAG_NR == 'NRP') then
                A_f = 0.0d0
            end if

            tStart3 = omp_get_wtime()
            call OMP_SET_NUM_THREADS(NTHREADS)
            !$OMP  parallel do&
            !$OMP& default (shared)&
            !$OMP& private (iel)
            do iel = 1, NEL_3d
                call FLOW_EQUATIONS(iel, FLAG_NR)
            end do
            !$OMP end parallel do
            tEnd = omp_get_wtime()
            if (output) write(*,'(a,t15,f9.4,a)',advance='no') 'Bulk time: ', tEnd-tStart3, ' s / '

            call DIRICHLET_BOUNDARY_CONDITIONS(1)
            call DIRICHLET_BOUNDARY_CONDITIONS(2)

            RES_NORM = dot_product(B_f,B_f)
            RES_NORM = sqrt(RES_NORM)

            tStart3 = omp_get_wtime()

            !LU DECOMPOSITION
            if (FLAG_NR == 'NRP') then

                PHASE_f = 22  ! ONLY FACTORIZATION

                call PARDISO&
                    (PT_f, MAXFCT_f, MNUM_f, MTYPE_f, PHASE_f, N_f, A_f, IA_f, CA_f,&
                    IDUM_f, NRHS_f, IPARM_f, MSGLVL_f, DDUM_f, DDUM_f, ERROR_f)

                if (ERROR_f /= 0) then
                    if (output) write(*,*)'THE FOLLOWING ERROR_F WAS DETECTED IN LU DECOMPOSITION: ', ERROR_f
                    stop
                end if

            end if

            !BACK SUBSTITUTION
            PHASE_f = 33    ! ONLY SOLUTION
            IPARM_f(8) = 2  ! MAX NUMBERS OF ITERATIVE REFINEMENT STEPS

            call PARDISO&
                (PT_f, MAXFCT_f, MNUM_f, MTYPE_f, PHASE_f, N_f, A_f, IA_f, CA_f,&
                IDUM_f, NRHS_f, IPARM_f, MSGLVL_f, B_f, S_f, ERROR_f)

            if (ERROR_f /= 0) then
                if (output) write(*,*)'THE FOLLOWING ERROR_F WAS DETECTED IN LU DECOMPOSITION: ', ERROR_f
                stop
            end if

            tEnd = omp_get_wtime()
            if (output) write(*,'(a,t15,f9.4,a)') 'LS time: ', tEnd-tStart3, ' s'

            call CHECK_CONVERGENCE(ITER_f, xNITER, FLAG_NR, COR_NORM_NEW_f, &
                COR_NORM_OLD_f, S_f, RES_NORM, LMSG, output)

            if (ITER_f == 1) Cor_Norm0 = COR_NORM_NEW_f

            if (LMSG) then
                EMERGENCY = .true.
                goto 300
            end if

            tStart3 = omp_get_wtime()

            k = 0
            do i = 1, NODTOL
                do j = 1, NEQ_f
                    k = k + 1
                    TL(i,j) = TL(i,j) - xF*S_f(k)
                end do
            end do

            ! call OMP_SET_NUM_THREADS(NTHREADS)
            !$OMP  parallel do&
            !$OMP& default (shared)&
            !$OMP& private (iel)
            do iel = 1, NEL_3d
                call STRESS_STORAGE( iel )
            end do
            !$OMP end parallel do

            tEnd = omp_get_wtime()

            if (output) write(*,'(a,t15,f9.4,a)') 'Store time: ', tEnd-tStart3, ' s'

            if (output) write(*,'(a,t15,f9.4,a)') 'NR time: ', tEnd-tStart2, ' s'
            if (output) write(*,'(a)') '----------------------------------------------------------------------'

        end do LOOP_loopNewtonRapshon

        tEnd = omp_get_wtime()
        if (output) write(*,'(a,t15,f9.4,a)') 'Solution time: ', tEnd-tStart, ' s'

    end subroutine loopNewtonRapshon

    !-----------------------------------------------------------------------

    subroutine CHECK_CONVERGENCE&
        (iter, maxiter, Flag_NR, Cor_Norm, Cor_norm_old, S_f, Res_norm, LMSG, output)

        ! use TIME_INTEGRATION_MODULE, only: Increment
        use NRAPSHON_MODULE, only: ERROR_NR, MNR_eps

        implicit none

        integer, intent(in) :: iter, maxiter
        character(*), intent(inout) :: Flag_NR
        real(8), intent(out) :: Cor_Norm
        real(8), intent(inout) :: Cor_norm_old
        real(8), dimension(:), intent(in) :: S_f
        real(8), intent(in) :: Res_norm
        logical, intent(out) :: LMSG
        logical, intent(in) :: output

        LMSG = .false.

        Cor_Norm  = sqrt(dot_product(S_f,S_f))

        if (Cor_Norm /= Cor_Norm) then
            if (output) write(*,'(a)') 'NaN ENCOUNTERED, GENERAL STOP'
            LMSG = .true.
            return
        end if

        KIND_OF_NEWTON_RAPHSON: SELECT case(Flag_NR)

        case ('NRP')

            if (output) write(*, 51)iter, Res_norm, Cor_Norm

            if ( Cor_Norm < MNR_eps ) then
                Flag_NR = 'MNR'
            else
                Flag_NR = 'NRP'
            end if

        case ('MNR')

            if (output) write(*, 52)iter, Res_norm, Cor_Norm

            if (Cor_Norm > Cor_norm_old) then
                Flag_NR = 'NRP'
            else
                Flag_NR = 'MNR'
            end if

            ! if (Increment <= 200) then
            !     if ( Cor_Norm  <=  ERROR_NR ) then
            !         Flag_NR = 'NRP'
            !     end if
            ! end if

        case default
            if (output) write(*,'(a)')' INCORRECT CHOICE OF NEWTON RAPHSON FLAG '
            stop
        end select KIND_OF_NEWTON_RAPHSON

        !  ONLY FULL NEWTON RAPHSON 
        ! Flag_NR = 'NRP'

        if ( Cor_Norm > 1.0d8 ) then
            if (output) write(*,'(a)')'TOO LARGE NORMA !!!'
            LMSG = .true.
            return
        end if


        if ( iter > maxiter ) then
            if (output) write(*,'(a)')'TOO MANY ITERATIONS !!!'
            LMSG = .true.
            return
        end if

        Cor_norm_old = Cor_Norm

        51 FORMAT(i5,'. NRP: ','Res Norm = ', es12.4, ',  Cor Norm = ', es12.4)
        52 FORMAT(i5,'. MNR: ','Res Norm = ', es12.4, ',  Cor Norm = ', es12.4)

    end subroutine CHECK_CONVERGENCE

end module Newton_mod
