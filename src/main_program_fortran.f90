program main_program_fortran
    use Quantum_Ising_1D
    use mkl_spblas
    implicit none

    integer(8) :: number_sites, nev, ncv
    real(8) :: trans_magnet_field, long_magnet_field, long_magnetiz, trans_magnetiz
    real(8), allocatable :: eval_arr(:), evec_matrix(:,:)

    ! real(8), allocatable :: in(:), out(:,:)
    ! integer :: i,stat

    number_sites = 19
    nev = 1
    ncv = 10

    long_magnet_field = 0
    trans_magnet_field = 1.1

    print *, 'IF YOU COMPILE WITH INTEL, TRUE=-1.'
    print *, 'THIS PROGRAM RELIES ON INTEL UGLY CONVENTION'
    print *, '  '

    call Initialize_Quantum_Ising(number_sites, nev, ncv)
    call Initialize_Hamiltonian(long_magnet_field, trans_magnet_field)
    ! allocate(in(tot_Hilbert_elem),out(tot_Hilbert_elem,tot_Hilbert_elem))

    ! do i=1,tot_Hilbert_elem
    !     in=0.d0
    !     in(i)=1.d0
    !     print*,i
    !     stat = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE,1.d0,sp_matrix,descr_sp_matrix,in,0.d0,out(i,:))
    !     print*, stat
    ! enddo
    ! call PRINT_MATRIX('',out)
    allocate(eval_arr(nev))
    allocate(evec_matrix(tot_Hilbert_elem, nev))

    call diag_lanczos('I', 'SA', eval_arr, evec_matrix)

    call get_long_magnetiz(tot_Hilbert_elem, evec_matrix(:,nev), long_magnetiz) !evec_matrix(nev-x) = evector of x state
    call get_trans_magnetiz(tot_Hilbert_elem, evec_matrix(:,nev), trans_magnetiz)

    print *, long_magnetiz
    print *, trans_magnetiz
    print *, eval_arr(nev)

    contains

    SUBROUTINE PRINT_MATRIX( DESC, A)
        CHARACTER*(*) :: DESC
        INTEGER :: M, N
        DOUBLE PRECISION ::  A(:,:)
        
        INTEGER :: I, J
        
        M=size(A,1)
        N=size(A,2)
        WRITE(100,*)
        WRITE(100,*) DESC
        DO I = 1, M
            WRITE(100,9998) ( A( I, J ), J = 1, N )
        END DO
        9998 FORMAT( 1000(:,1X,F10.6) )
        RETURN
        END


end program main_program_fortran