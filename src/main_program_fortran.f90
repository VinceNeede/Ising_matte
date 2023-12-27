program main_program_fortran
    use Quantum_Ising_1D
    implicit none

    integer(8) :: number_sites, nev, ncv
    real(8) :: trans_magnet_field, long_magnet_field, long_magnetiz, trans_magnetiz
    real(8), allocatable :: eval_arr(:), evec_matrix(:,:)

    number_sites = 22
    nev = 1
    ncv = 10

    long_magnet_field = 0.d0
    trans_magnet_field = 1.1d0

    print *, 'IF YOU COMPILE WITH INTEL, TRUE=-1.'
    print *, 'THIS PROGRAM RELIES ON INTEL UGLY CONVENTION'
    print *, '  '

    call Initialize_Quantum_Ising(number_sites, nev, ncv)
    call Initialize_Hamiltonian(long_magnet_field, trans_magnet_field)

    allocate(eval_arr(nev))
    allocate(evec_matrix(tot_Hilbert_elem, nev))

    call diag_lanczos('I', 'SA', eval_arr, evec_matrix)

    call get_long_magnetiz(tot_Hilbert_elem, evec_matrix(:,nev), long_magnetiz) !evec_matrix(nev-x) = evector of x state
    call get_trans_magnetiz(tot_Hilbert_elem, evec_matrix(:,nev), trans_magnetiz)

    print *, long_magnetiz
    print *, trans_magnetiz
    print *, eval_arr(nev)


end program main_program_fortran