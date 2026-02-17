
module Tools_mod

    implicit none

    contains

    function traceTensor(A)

        implicit none

        real(8), dimension(:,:), intent(in) :: A
        real(8) :: traceTensor

        integer :: i

        traceTensor = 0.0d0
        do i = 1, size(A,1)
            traceTensor = traceTensor + A(i,i)
        end do

    end function traceTensor

    !-----------------------------------------------------------------------  

    function doubleDotProduct(A, B)

        implicit none

        real(8), dimension(:,:), intent(in) :: A, B
        real(8) :: doubleDotProduct

        integer :: i, j

        doubleDotProduct = 0.0d0
        do i = 1, size(A,1)
            do j = 1, size(A,2)
                doubleDotProduct = doubleDotProduct + A(i,j)*B(j,i)
            end do
        end do

    end function doubleDotProduct

    !----------------------------------------------------------------------- 

    function tensorMagnitude(A)

        implicit none

        real(8), dimension(:,:), intent(in) :: A
        real(8) :: TensorMagnitude

        real(8) :: double_dot_A

        double_dot_A = doubleDotProduct(A,A)

        TensorMagnitude = sqrt(0.5d0*double_dot_A)

    end function tensorMagnitude

    !-----------------------------------------------------------------------  

    function determinant2D(A)

        implicit none

        real(8), dimension(2,2), intent(in) :: A
        real(8) :: determinant2D

        determinant2D = A(1,1)*A(2,2) - A(1,2)*A(2,1)

    end function determinant2D

    !--------------------------------------------------------------------------

    function determinant3D(A)

        implicit none

        real(8), dimension(3,3), intent(in) :: A
        real(8) :: determinant3D

        real(8), dimension(2,2) :: temp_A
        real(8) :: term1, term2, term3

        temp_A(:,:) = A(2:3,2:3)
        term1 = A(1,1)*determinant2D(temp_A)

        temp_A(:,:) = A(2:3,1:3:2)
        term2 = A(1,2)*determinant2D(temp_A)

        temp_A(:,:) = A(2:3,1:2)
        term3 = A(1,3)*determinant2D(temp_A)

        determinant3D = term1-term2+term3

    end function determinant3D
    
    !-----------------------------------------------------------------------  

    subroutine inverseTensor1D(A, AI)

        implicit none

        real(8), dimension(1,1), intent(in) :: A
        real(8), dimension(1,1), intent(out) :: AI

        AI(1,1) =  1.0d0/A(1,1)
        
    end subroutine inverseTensor1D

    !-----------------------------------------------------------------------  

    subroutine inverseTensor2D(A, AI)

        implicit none

        real(8), dimension(:,:), intent(in) :: A
        real(8), dimension(:,:), intent(out) :: AI

        real(8) :: detA

        detA = determinant2D(A)

        AI(:,:) = 0.0d0
        AI(1,1) =  A(2,2) ; AI(1,2) = -A(1,2)
        AI(2,1) = -A(2,1) ; AI(2,2) =  A(1,1)

        AI(:,:) = AI(:,:)/detA

    end subroutine inverseTensor2D

    !-----------------------------------------------------------------------  

    subroutine inverseTensor3D(A, AI)

        implicit none

        real(8), dimension(:,:), intent(in) :: A
        real(8), dimension(:,:), intent(out) :: AI

        ! real(8) :: detA
        ! real(8), dimension(2,2) :: temp_A
        ! integer :: i, j
        ! integer, dimension(3,3), parameter :: &
        !     signs = reshape( [+1,-1,+1,-1,+1,-1,+1,-1,+1], shape(signs) )
        ! integer, dimension(2,3), parameter :: &
        !     irows = reshape( [2,3, 1,3, 1,2], shape = shape(irows) )
        real(8) :: detS, det11, det12, det13, det21, det22, det23, det31, det32, det33

        ! do i = 1, 3
        !     do j = 1, 3
        !         temp_A(:,:) = A(irows(:,i),irows(:,j))
        !         AI(i,j) = signs(i,j)*determinant2D(temp_A)
        !     end do
        ! end do

        ! AI(:,:) = transpose(AI)

        ! detA = determinant3D(A)

        ! AI(:,:) = AI(:,:)/detA


        det11 = A(2,2)*A(3,3) - A(2,3)*A(3,2)
        det12 = A(2,1)*A(3,3) - A(2,3)*A(3,1)
        det13 = A(2,1)*A(3,2) - A(2,2)*A(3,1)

        det21 = A(1,2)*A(3,3) - A(1,3)*A(3,2)
        det22 = A(1,1)*A(3,3) - A(1,3)*A(3,1)
        det23 = A(1,1)*A(3,2) - A(1,2)*A(3,1)

        det31 = A(1,2)*A(2,3) - A(1,3)*A(2,2)
        det32 = A(1,1)*A(2,3) - A(1,3)*A(2,1)
        det33 = A(1,1)*A(2,2) - A(1,2)*A(2,1)

        detS = A(1,1)*det11 - A(1,2)*det12 + A(1,3)*det13

        AI(:,:) = 0.0d0
        AI(1,1) =  det11 ; AI(1,2) = -det12 ; AI(1,3) =  det13
        AI(2,1) = -det21 ; AI(2,2) =  det22 ; AI(2,3) = -det23
        AI(3,1) =  det31 ; AI(3,2) = -det32 ; AI(3,3) =  det33

        AI = transpose(AI)/detS

    end subroutine inverseTensor3D
    
    !---------------------------------------------------------------------

    subroutine getStressComponent(i, istart, j, k)
        
        implicit none

        integer, intent(in) :: i, istart
        integer, intent(out) :: j, k

        if (i == istart) then
            j = 1 ; k = 1
        else if (i == istart+1) then
            j = 1 ; k = 2
        else if (i == istart+2) then
            j = 1 ; k = 3
        else if (i == istart+3) then
            j = 2 ; k = 2
        else if (i == istart+4) then
            j = 2 ; k = 3
        else if (i == istart+5) then
            j = 3 ; k = 3
        else
            write(*,'(a)') 'Out of range getStressComponent!'
            stop
        end if

    end subroutine getStressComponent
    
    !---------------------------------------------------------------------

    integer function getParameterID(par) result(id)

        implicit none

        character(*), intent(in) :: par

        select case (par)

        case ('TdN')
            id = 1
        case ('GveN')
            id = 2
        case ('GeN')
            id = 3
        case ('JmN', 'BetaN')
            id = 4
        case ('EpN')
            id = 5
        case ('BetaN2', 'TaN')
            id = 6
        case default
            write(*,'(a)') 'Wrong parameter name in getParameterID!'
            stop
        end select

    end function getParameterID

    !----------------------------------------------------------------------

    subroutine basis3D(Xi_el, bfn, dfdci, dfdxi, Xi_gs, dXidci, CJAC)

        implicit none

        real(8), dimension(:,:), intent(in) :: Xi_el, dfdci
        real(8), dimension(:), intent(in) :: bfn
        real(8), dimension(:,:), intent(out) :: dfdxi, dXidci
        real(8), dimension(:), intent(out) :: Xi_gs
        real(8), intent(out) :: CJAC

        integer :: i, Nbf
        real(8) :: Cx, Cy, Cz, Ex, Ey, Ez, Sx, Sy, Sz

        Nbf = size(Xi_el,1)
        Xi_gs(:) = 0.0d0 ; dXidci(:,:) = 0.0d0
        do i = 1, Nbf
            Xi_gs(1) = Xi_gs(1) + bfn(i)*Xi_el(i,1)
            dXidci(1,1) = dXidci(1,1) + dfdci(i,1)*Xi_el(i,1)
            dXidci(1,2) = dXidci(1,2) + dfdci(i,2)*Xi_el(i,1)
            dXidci(1,3) = dXidci(1,3) + dfdci(i,3)*Xi_el(i,1)

            Xi_gs(2) = Xi_gs(2) + bfn(i)*Xi_el(i,2)
            dXidci(2,1) = dXidci(2,1) + dfdci(i,1)*Xi_el(i,2)
            dXidci(2,2) = dXidci(2,2) + dfdci(i,2)*Xi_el(i,2)
            dXidci(2,3) = dXidci(2,3) + dfdci(i,3)*Xi_el(i,2)

            Xi_gs(3) = Xi_gs(3) + bfn(i)*Xi_el(i,3)
            dXidci(3,1) = dXidci(3,1) + dfdci(i,1)*Xi_el(i,3)
            dXidci(3,2) = dXidci(3,2) + dfdci(i,2)*Xi_el(i,3)
            dXidci(3,3) = dXidci(3,3) + dfdci(i,3)*Xi_el(i,3)
        end do

        associate(dXdC => dXidci(1,1), dXdE => dXidci(1,2), dXdS => dXidci(1,3), &
            dYdC => dXidci(2,1), dYdE => dXidci(2,2), dYdS => dXidci(2,3), &
            dZdC => dXidci(3,1), dZdE => dXidci(3,2), dZdS => dXidci(3,3))

            CJAC = dXdC*(dYdE*dZdS-dYdS*dZdE) - &
                    dXdE*(dYdC*dZdS-dYdS*dZdC) + &
                    dXdS*(dYdC*dZdE-dYdE*dZdC)

            Cx = (dYdE*dZdS-dYdS*dZdE)/CJAC
            Cy =-(dXdE*dZdS-dZdE*dXdS)/CJAC
            Cz = (dXdE*dYdS-dYdE*dXdS)/CJAC

            Ex =-(dYdC*dZdS-dZdC*dYdS)/CJAC
            Ey = (dXdC*dZdS-dZdC*dXdS)/CJAC
            Ez =-(dXdC*dYdS-dYdC*dXdS)/CJAC

            Sx = (dYdC*dZdE-dZdC*dYdE)/CJAC
            Sy =-(dXdC*dZdE-dZdC*dXdE)/CJAC
            Sz = (dXdC*dYdE-dYdC*dXdE)/CJAC

            do i = 1, Nbf
                dfdxi(i,1) = dfdci(i,1)*Cx + dfdci(i,2)*Ex + dfdci(i,3)*Sx 
                dfdxi(i,2) = dfdci(i,1)*Cy + dfdci(i,2)*Ey + dfdci(i,3)*Sy 
                dfdxi(i,3) = dfdci(i,1)*Cz + dfdci(i,2)*Ez + dfdci(i,3)*Sz 
            end do
        end associate

    end subroutine basis3D

    !---------------------------------------------------------------------

    subroutine writeElapsedTime(start_time,end_time,str)

        implicit none

        real(8), intent(in) :: start_time, end_time
        character(len=*), intent(in) :: str

        real(8) :: total_time, seconds
        integer  :: days, hours, minutes

        total_time = end_time-start_time
        days    = int(total_time/8.64d4)
        hours   = int((total_time-days*8.64d4)/3.6d3)
        minutes = int((total_time-days*8.64d4-hours*3.6d3)/60.0d0)
        seconds = total_time-days*8.64d4-hours*3.6d3-minutes*60.0d0

        write(*, 60) str//': ', days, hours, minutes, seconds
        write(*, '(a)') '----------------------------------------------------------------------'

        60 format(A, T30, I3, " d, ", I3, " hr, ", I3, " min, ", F6.3, " s")

    end subroutine writeElapsedTime
    
    !-----------------------------------------------------------------------

    subroutine extrapolationLagrange(Xp, Xb, Xo, X, Lb, Lo, L)

        implicit none

        real(8), intent(in) :: Xp, Xb, Xo, X
        real(8), intent(out) :: Lb, Lo, L

        Lb = (Xp-Xo)*(Xp-X) /((Xb-Xo)*(Xb-X))
        Lo = (Xp-Xb)*(Xp-X) /((Xo-Xb)*(Xo-X))
        L  = (Xp-Xb)*(Xp-Xo)/((X -Xb)*(X -Xo))

    end subroutine extrapolationLagrange

    !-----------------------------------------------------------------------

    real(8) function perturbVariable ( x ) result(eps)

        implicit none

        real(8), intent(in) :: x

        real(8) :: EP_RES = 1.0d-9

        eps = EP_RES*DMAX1(1.0d0,abs(x))

    end function perturbVariable

    !----------------------------------------------------------------------

    subroutine storeResidual( TEMP, NM, IDM, INEQ, B, IDIM_B )

        implicit none

        integer, intent(in) :: IDM, INEQ, IDIM_B

        integer, dimension(IDM), intent(in) :: NM
        real(8), dimension(IDM,INEQ), intent(in) :: TEMP
        real(8), dimension(IDIM_B), intent(inout) :: B

        integer :: I, IEQ, IROW

        do I = 1, IDM
            do IEQ = 1, INEQ
                IROW  = NM(I) + IEQ - 1
                !$OMP ATOMIC
                B(IROW) = B(IROW) + TEMP(I,IEQ)
            end do
        end do

    end subroutine storeResidual

    !----------------------------------------------------------------------

    subroutine storeJacobian(TP, IDM, JDIM, JNEQ, INEQ, NM, IA, IIA, CSR, ICSR, A, IDIM_A)

        implicit none

        integer, intent(in) :: ICSR, IIA
        integer, intent(in) :: IDM, JDIM, INEQ, JNEQ
        integer, intent(in) :: IDIM_A

        integer, dimension(ICSR), intent(in) :: CSR
        integer, dimension(IIA),  intent(in) :: IA
        integer, dimension(IDM),  intent(in) :: NM

        real(8), dimension(IDM,JDIM,JNEQ,INEQ), intent(in)    :: TP
        real(8), dimension(IDIM_A),             intent(inout) :: A

        integer :: I, J, L, IEQ, JEQ
        integer :: IROW, ICOL, IAD

        do I = 1, IDM
            IROW = NM(I)
            do IEQ = 1, INEQ
                IAD  = (IEQ-1)*(IA(IROW+1)-IA(IROW))
                do J = 1, JDIM
                    L = (I-1)*JDIM + J
                    do JEQ = 1, JNEQ
                    ICOL = IAD + CSR(L) + JEQ - 1
                        !$OMP ATOMIC
                        A(ICOL) = A(ICOL) + TP(I,J,JEQ,IEQ)
                    end do
                end do
            end do
        end do

    end subroutine storeJacobian

    !------------------------------------------------------------------------------

    subroutine readData(filename, x, y, z)

        implicit none

        character(len=*), intent(in)  :: filename
        real(8), allocatable, dimension(:), intent(out) :: x, y, z

        integer :: k, ierror, i, u
        real :: x_dump, y_dump, z_dump

        u = openReadFile(filename)

        k = 0
        do
            read(u,*,iostat=ierror) x_dump, y_dump, z_dump
            if (ierror /= 0) exit
            k = k+1
        end do

        
        if (allocated(x)) deallocate(x) ; allocate(x(k))
        if (allocated(y)) deallocate(y) ; allocate(y(k))
        if (allocated(z)) deallocate(z) ; allocate(z(k))

        rewind(u)

        do i = 1, k
            read(u,*,iostat=ierror) x(i), y(i), z(i)
        end do

        close(u)

    end subroutine readData

    !------------------------------------------------------------------------------

    integer function openReadFile(filename)

        implicit none

        character(len=*), intent(in) :: filename

        integer :: ierror, u

        u = 8
        open(unit=u, file = filename, action = 'read', iostat = ierror, position = 'rewind')
        if (ierror /= 0) then
            write(*,'(2A)') 'Cannot open ', filename
            stop
        end if

        openReadFile = u

    end function openReadFile

    !------------------------------------------------------------------------------

    logical function plotGraphs(iter)

        implicit none

        integer, intent(in) :: iter

        character(50) :: str

        ! call execute_command_line ("python3 plot.py")
        call execute_command_line ("python3.8 plot.py")
        write(str,'(i5)') iter
        str = 'cp -r Graphs out/Graphs_full/Graphs_'//trim(adjustl(str))
        call execute_command_line (str)

        plotGraphs = .true.

    end function plotGraphs

    !-----------------------------------------------------------------------

    subroutine DECOMP( NDIM, N, A, COND, IPVT, WORK )

        implicit none

        integer, intent(in) :: NDIM, N
        real(8), intent(inout) :: COND
        real(8), dimension(NDIM,N), intent(out) :: A
        real(8), dimension(N), intent(inout) :: WORK
        integer, dimension(N), intent(out) :: IPVT

        real(8) :: EK, T, ANORM, YNORM, ZNORM
        integer :: NM1, I, J, K, KP1, KB, KM1, M

        !    DECOMPOSES A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION
        !    AND ESTIMATES THE CONDITION OF THE MATRIX.
        !
        !    use SOLVE TO COMPUTE SOLUTIONS TO LINEAR SYSTEMS.
        !
        !    INPUT..
        !
        !    NDIM = DECLARED ROW DIMENSION OF THE ARRAY CONTAINING  A.
        !
        !    N = ORDER OF THE MATRIX.
        !
        !    A = MATRIX TO BE TRIANGULARIZED.
        !
        !    OUTPUT..
        !
        !    A   CONTAINS AN UPPER TRIANGULAR MATRIX  U  AND A PERMUTED
        !        VERSION OF A LOWER TRIANGULAR MATRIX  I-L  SO THAT
        !        (PERMUTATION MATRIX)*A = L*U .
        !
        !    COND = AN ESTIMATE OF THE CONDITION OF  A .
        !           FOR THE LINEAR SYSTEM  A*X = B, CHANGES IN  A  AND  B
        !           MAY CAUSE CHANGES  COND  TIMES AS LARGE IN  X .
        !           if  COND+1.0 == COND , A IS SINGULAR TO WORKING
        !           PRECISION.  COND  IS SET TO  1.0D+32  if EXACT
        !           SINGULARITY IS DETECTED.
        !
        !    IPVT = THE PIVOT VECTOR.
        !           IPVT(K) = THE INDEX OF THE K-TH PIVOT ROW
        !           IPVT(N) = (-1)**(NUMBER OF INTERCHANGES)
        !
        !    WORK SPACE..  THE VECTOR  WORK  MUST BE DECLARED AND INCLUDED
        !           IN THE CALL.  ITS INPUT CONTENTS ARE IGNORED.
        !           ITS OUTPUT CONTENTS ARE USUALLY UNIMPORTANT.
        !
        !    THE DETERMINANT OF A CAN BE OBTAINED ON OUTPUT BY
        !           DET(A) = IPVT(N) * A(1,1) * A(2,2) * ... * A(N,N).
        !

        IPVT(N) = 1
        if (N == 1) GO TO 80
        NM1 = N - 1

        !    COMPUTE 1-NORM OF A
        ANORM = 0.0D0

        do J = 1, N
            T = 0.0D0
            do I = 1, N
                T = T + abs(A(I,J))
            end do
            if (T > ANORM) ANORM = T
        end DO

        !    GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
        LOOP_35: do K = 1,NM1

            KP1 = K + 1

            !      FIND PIVOT
            M = K
            do I = KP1, N
                if (abs(A(I,K)) > abs(A(M,K))) M = I
            end do

            IPVT(K) = M
            if (M /= K) IPVT(N) = -IPVT(N)
            T = A(M,K)
            A(M,K) = A(K,K)
            A(K,K) = T

            !      SKIP STEP if PIVOT IS ZERO
            if (T == 0.0D0) cycle LOOP_35

            !      COMPUTE MULTIPLIERS
            do I = KP1, N
                A(I,K) = -A(I,K) / T
            end do

            !      INTERCHANGE AND ELIMINATE BY COLUMNS
            LOOP_30: do J = KP1,N
                T = A(M,J)
                A(M,J) = A(K,J)
                A(K,J) = T
                if (T == 0.0D0) cycle LOOP_30
                do I = KP1, N
                    A(I,J) = A(I,J) + A(I,K)*T
                end do
            end do LOOP_30

        end do LOOP_35

        !    COND = (1-NORM OF A)*(AN ESTIMATE OF 1-NORM OF A-INVERSE)
        !    ESTIMATE OBTAINED BY ONE STEP OF INVERSE ITERATION FOR THE
        !    SMALL SINGULAR VECTOR.  THIS INVOLVES SOLVING TWO SYSTEMS
        !    OF EQUATIONS, (A-TRANSPOSE)*Y = E  AND  A*Z = Y  WHERE  E
        !    IS A VECTOR OF +1 OR -1 CHOSEN TO CAUSE GROWTH IN Y.
        !    ESTIMATE = (1-NORM OF Z)/(1-NORM OF Y)
        !
        !    SOLVE (A-TRANSPOSE)*Y = E

        do K = 1, N
            T = 0.0D0
            if (K == 1) GO TO 45
            KM1 = K-1
            do I = 1, KM1
                T = T + A(I,K)*WORK(I)
            end do
            45 EK = 1.0D0
            if (T < 0.0D0) EK = -1.0D0
            if (A(K,K) == 0.0D0) GO TO 90
            WORK(K) = -(EK + T)/A(K,K)
        end do


        LOOP_60: do KB = 1, NM1
            K = N - KB
            T = 0.0D0
            KP1 = K+1
            do I = KP1, N
                T = T + A(I,K)*WORK(I)
            end do
            WORK(K) = T + WORK(K)
            M = IPVT(K)
            if (M == K) cycle LOOP_60
            T = WORK(M)
            WORK(M) = WORK(K)
            WORK(K) = T
        end do LOOP_60

        !
        YNORM = 0.0D0
        do I = 1, N
            YNORM = YNORM + abs(WORK(I))
        end do

        !    SOLVE A*Z = Y
        call SOLVE( NDIM, N, A, WORK, IPVT )

        ZNORM = 0.0D0
        do I = 1, N
            ZNORM = ZNORM + abs(WORK(I))
        end do

        !    ESTIMATE CONDITION
        COND = ANORM*ZNORM/YNORM
        if (COND < 1.0D0) COND = 1.0D0
        return

        !    1-BY-1
        80 COND = 1.0D0
        if (A(1,1) /= 0.0D0) return

        !    EXACT SINGULARITY
        90 COND = 1.0D+32

    end subroutine DECOMP

    !-----------------------------------------------------------------------

    subroutine SOLVE( NDIM, N, A, B, IPVT )

        implicit none

        integer, intent(in) :: NDIM, N
        integer, dimension(N), intent(in)  :: IPVT
        real(8), dimension(N), intent(out) :: B
        real(8), dimension(NDIM,N), intent(in)  :: A

        integer :: KB, KM1, NM1, KP1, I, K, M
        real(8) :: T

        !    SOLUTION OF LINEAR SYSTEM, A*X = B .
        !    do NOT use if DECOMP HAS DETECTED SINGULARITY.
        !
        !    INPUT..
        !
        !    NDIM = DECLARED ROW DIMENSION OF ARRAY CONTAINING A .
        !
        !    N = ORDER OF MATRIX.
        !
        !    A = TRIANGULARIZED MATRIX OBTAINED FROM DECOMP .
        !
        !    B = RIGHT HAND SIDE VECTOR.
        !
        !    IPVT = PIVOT VECTOR OBTAINED FROM DECOMP .
        !
        !    OUTPUT..
        !
        !    B = SOLUTION VECTOR, X .

        !    FORWARD ELIMINATION
        if (N == 1) go to 50

        NM1 = N-1

        do K = 1, NM1
            KP1 = K+1
            M = IPVT(K)
            T = B(M)
            B(M) = B(K)
            B(K) = T
            do I = KP1, N
                B(I) = B(I) + A(I,K)*T
            end do
        end do

        !    BACK SUBSTITUTION
        do KB = 1,NM1
            KM1 = N-KB
            K = KM1+1
            B(K) = B(K)/A(K,K)
            T = -B(K)
            do I = 1, KM1
                B(I) = B(I) + A(I,K)*T
            end do
        end do

        50 B(1) = B(1)/A(1,1)

    end subroutine SOLVE

end module Tools_mod
