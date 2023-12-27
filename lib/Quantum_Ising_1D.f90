module Quantum_Ising_1D
    use mkl_spblas
    implicit none
    private
    public Initialize_Quantum_Ising, Initialize_Hamiltonian, get_long_magnetiz, get_trans_magnetiz, diag_lanczos, tot_Hilbert_elem

    integer, parameter :: number_spin_states = 2
    integer, parameter :: dimensions = 1
    integer, parameter :: boundary_conditions_flag = 1
    integer :: number_sites, total_sites, tot_Hilbert_elem, nev, ncv, ldv, nnz
    logical, allocatable :: comput_base(:,:)
    real(8),allocatable :: val_arr(:)
    integer(8),allocatable :: col_arr(:), row_arr(:)

    type(sparse_matrix_T) :: sp_matrix
    type(MATRIX_DESCR) :: descr_sp_matrix !caratteristiche della matrice, come sym, hermitiana, ecc

    ! FOR ALL THE FOLLOWING LOCAL NAMES, THEY ARE LEFT AS IN THE ORIGINAL DOCUMENTATION
    ! Allocations are done in initialize_comput_base
 
    !     %------------------------------------------------------%
    !     | NCV is the largest number of basis vectors that will |
    !     |     be used in the Implicitly Restarted Arnoldi      |
    !     |     Process.  Work per major iteration is            |
    !     |     proportional to N*NCV*NCV.                       |
    !     |                                                      |
    !     |NEV is the number of eigenvalues                      |
    !     %------------------------------------------------------%

    !     %------------------------------------------------------%
    !     |  Note: NEV and NCV must satisfy the following        |
    !     | conditions:                                          |
    !     |                  NEV + 1 <= NCV                      |
    !     %------------------------------------------------------%

    !     %--------------%
    !     | Local Arrays |
    !     %--------------%
   
    Real(8), allocatable :: lancz_basis_matrix(:,:), workl(:), workd(:), resid(:), ax(:)
    logical, allocatable :: select(:)
    integer :: iparam(11), ipntr(11)
    
    !     %---------------%
    !     | Local Scalars |
    !     %---------------%
    
    integer(4) :: previous_file_digits
    integer :: ido, length_workl, info, ierr, j, ishfts, maxitr, mode1, nconv, itry
    Real(8) :: tol, sigma, rnorm    !tol(arpack.ng) == tol(Spectra C++), 
                                    !maxitr(arpack.ng) == maxit(Spectra C++)
    
    !     %------------%
    !     | Parameters |
    !     %------------%
    
    Real(8), parameter :: zero=0.0d+0  !machine precision, will be used for tol
    
    !     %-----------------------------%
    !     | BLAS & LAPACK routines used |
    !     %-----------------------------%   

    !   %------------------------------------------------------------%
    !   |                       Routines called:                     |
    !   | dsaupd  ARPACK reverse communication interface routine.    |
    !   | dseupd  ARPACK routine that returns Ritz values            |
    !   |       and (optionally) Ritz vectors.                       |
    !   | dnrm2   Level 1 BLAS that computes the norm of a vector.   |
    !   | daxpy   Level 1 BLAS that computes y <- alpha*x+y.         |
    !   %------------------------------------------------------------%
    Real(8) dnrm2  
    external dnrm2, daxpy, dseupd, dsaupd, dcopy, dgetv0


    


    contains

    subroutine Initialize_Quantum_Ising(number_sites_external, nev_external, ncv_external)
        implicit none
        integer, intent(in) :: number_sites_external, nev_external, ncv_external
        integer :: i_row, j_col, i_holder ! j_holder

        number_sites = number_sites_external
        nev = nev_external
        ncv = ncv_external
        
        total_sites = number_sites*dimensions
        tot_Hilbert_elem = number_spin_states**total_sites
        ldv = tot_Hilbert_elem
        nnz = tot_Hilbert_elem*(total_sites+1)

        !Inizializing all indipendent configuration of Hilbert Space.

        ! %-----------------------------------------------------%
        ! |             OPTIMIZATION TURNED OFF                 |    
        ! |                  (TRANSPOSITION)                    |
        ! | column major language, comput_base has consecutive  |
        ! | calls on adjacent sites, not hilbert states         |
        ! %-----------------------------------------------------%
       
        ! allocate(comput_base(total_sites, tot_Hilbert_elem)) !! TRASPOSE FROM C++
        ! j_holder = 0
        ! do j_col = 1, tot_Hilbert_elem !column major language
        !     j_holder = j_col
        !     do i_row = 1, total_sites
        !         comput_base(i_row, j_col) = MOD(j_holder,2)
        !         j_holder = j_holder/2
        !     end do
        ! end do

        allocate(comput_base(tot_Hilbert_elem, total_sites)) !! SAME AS C++
        i_holder = 0
        do i_row = 1, tot_Hilbert_elem
            i_holder = i_row - 1
            do j_col = 1, total_sites
                comput_base(i_row, j_col) = MOD(i_holder,2)
                i_holder = i_holder/2
            end do
        end do

    !     %-------------------------%
    !     | Local Arrays allocation |
    !     %-------------------------%

        allocate(lancz_basis_matrix(ldv,ncv)) 
        allocate(row_arr(nnz),col_arr(nnz),val_arr(nnz))
        length_workl = ncv*(ncv+8)
        allocate(workl(length_workl))
        allocate(workd(3*tot_Hilbert_elem))
        allocate(resid(tot_Hilbert_elem))
        allocate(ax(tot_Hilbert_elem))
        allocate(select(ncv)) 

    end subroutine Initialize_Quantum_Ising


    subroutine Initialize_Hamiltonian(long_magnet_field, trans_magnet_field)
        implicit none
        real(8), intent(in) :: long_magnet_field, trans_magnet_field
        integer :: i_row, j_col !indices for computational base we are scrolling through.
        integer :: new_i_row    !index for non diagonal elements of sp_matrix.
        integer :: index_holder
        integer :: stat, try           
        stat = 0                !If altered an error occurred

        val_arr = 0.d0
        
        descr_sp_matrix % TYPE = sparse_matrix_TYPE_SYMMETRIC            !!! col .geq. row
        descr_sp_matrix % Mode = SPARSE_FILL_MODE_UPPER
        descr_sp_matrix % diag = SPARSE_DIAG_NON_UNIT

    ! %-------------------------------------------------------------------------------------%
    ! |  we are finding the hamiltonian value for each configuration (cfr comput_base).     |
    ! |  Each one needs to valuate Spin values over each of the sites.                      |
    ! |                                                                                     |
    ! |  We build a matrix with, as eigenvalues (i.e. diagonals elements if diagonalized),  |
    ! |  all these hamiltonians.                                                            |
    ! |                                                                                     |
    ! |  In conclusion we find the configuration which_evalmizes the energy,                |                                          
    ! |  associated with the hamiltonian which_evalmizes the enrgy itself.                  |
    ! |  However trasverse field introduces non diagonal elements in the sp_matrix          |  
    ! %-------------------------------------------------------------------------------------%

    ! %---------------------------------------------------------------------------------%
    ! |  We place coupling term along with longitudinal.                                |
    ! |  Defining Pauli matrixes as Sigma_z=(1,0,-1,0),                                 |
    ! |  then we shall put the longitudinal directly on z,                              |
    ! |  so that we can approach to an already diagonalized term.                       |
    ! |                                                                                 |
    ! |  Sigma_z adds to the Hamiltonian  a value that depends on the spin of the site  |
    ! |  " 2 different values for Sigma_z(0,0) and Sigma_z(1,1) "                       |
    ! |                                                                                 |
    ! |  Sigma_x and Sigma_y adds a constant to the Hamiltonian                         |
    ! |  and flips the state in the Hilbert space.                                      | 
    ! |  Doing so is equal to change the configuration itself.                          |
    ! %---------------------------------------------------------------------------------%
    
    ! %---------------------------------------------------------%
    ! |  matrix:  0 1                                           |
    ! |           2 3                                           |
    ! |                                                         |
    ! | row_arr = (1, 2, 2)                                     |
    ! | col_arr = (2, 1, 2)                                     |
    ! | val_arr = (1, 2, 3)                                     |
    ! |                                                         |
    ! | or:                                                     |
    ! |                                                         |
    ! | row_arr = (2, 2, 1)                                     |
    ! | col_arr = (2, 1, 2)                                     |
    ! | val_arr = (3, 2, 1)                                     |
    ! |                                                         |
    ! | val = 0 are not considered                              |
    ! |                                                         |
    ! | row_arr(i), col_arr(i) are row and column of val(i)     |
    ! | It's not important in what order you put values in arr  |
    ! | It's important that (i) is the same for all 3 arrays    |
    ! |                                                         |
    ! |                                                         |
    ! | The following array order will be:                      |
    ! | diagonal elements as first tot_Hilbert_elem             |
    ! | non diagonal as the other ones                          |
    ! %---------------------------------------------------------%

        ! fortran is a column major language

        !   %---------------------------------------%
        !   |   (Sigma_(z) * Sigma_(z+1)) coupling  |
        !   %---------------------------------------%

        do i_row=1, tot_Hilbert_elem
            col_arr(i_row) = i_row
            row_arr(i_row) = i_row
            do j_col=1, total_sites-1

                if (comput_base(i_row, j_col) == comput_base(i_row, j_col+1)) then
                    val_arr(i_row) = val_arr(i_row) - 1.d0
                else
                    val_arr(i_row) = val_arr(i_row) + 1.d0
                end if

            end do
        end do

        if (boundary_conditions_flag==1) then
            do i_row=1, tot_Hilbert_elem
                ! col_arr(i_row) = i_row
                ! row_arr(i_row) = i_row

                if (comput_base(i_row, total_sites) == comput_base(i_row, 1)) then
                    val_arr(i_row) = val_arr(i_row) - 1.d0
                else
                    val_arr(i_row) = val_arr(i_row) + 1.d0
                end if

            end do
        end if

        
        !   %-----------------------%
        !   |   Longitudinal Field  |
        !   %-----------------------%

        do i_row=1, tot_Hilbert_elem
            ! col_arr(i_row) = i_row
            ! row_arr(i_row) = i_row

            do j_col=1, total_sites
                val_arr(i_row) = val_arr(i_row) + (1+2*comput_base(i_row, j_col)) * long_magnet_field
            end do
            ! if comput_base(i_row, j_col) == 1 then
            !     sp_matrix(i_row, j) = sp_matrix(i_row, j_col) - 1
            ! else 
            !     sp_matrix(i_row, j_col) = sp_matrix(i_row, j_col) + 1
            ! end if
        end do

        index_holder =  tot_Hilbert_elem !we have filled all diagonal elements 


        !   %---------------------%
        !   |   Transverse Field  |
        !   %---------------------%

        do i_row=1, tot_Hilbert_elem
            do j_col=1, total_sites
                index_holder = index_holder + 1
                new_i_row = i_row + (1+2*comput_base(i_row, j_col)) * 2**(j_col-1)

                col_arr(index_holder) = i_row
                row_arr(index_holder) = new_i_row
                val_arr(index_holder) = - trans_magnet_field ! val_arr(index_holder) - trans_magnet_field they do not overlap
            end do
            ! if comput_base(i_row, j_col) == 1 then
            !     new_i_row = i_row - 2**(j_col-1)
            ! else 
            !     new_i_row = i_row + 2**(j_col-1)
            ! end if
        end do     

        ! I tried to print all the values, they work.

        ! do i_row=1, tot_Hilbert_elem
        !     do j_col=1, total_sites
        !         print *, comput_base(i_row, j_col)
        !     end do
        !     print *, ' '
        ! end do 

        ! do try=1, nnz
        !     print *, int(val_arr(try)), row_arr(try), col_arr(try)
        ! end do
        ! print *, ''

        stat = mkl_sparse_d_create_coo(sp_matrix, SPARSE_INDEX_BASE_ONE, tot_Hilbert_elem, tot_Hilbert_elem,&
        nnz, row_arr, col_arr, val_arr)
        if (stat .ne. 0) print*, 'error in sparse matrix creation: ', stat   
    
    end subroutine Initialize_Hamiltonian



    subroutine get_long_magnetiz(len_evector, evector, long_magnetiz)
        implicit none
        integer, intent(in) :: len_evector
        real(8), intent(in) :: evector (len_evector)   
        real(8), intent(out) :: long_magnetiz
        
        integer :: i_row, j_col
        real(8):: sum_expect_val_sigma_z

        do i_row=1, tot_Hilbert_elem
            sum_expect_val_sigma_z = 0
            do j_col=1, total_sites
                sum_expect_val_sigma_z = sum_expect_val_sigma_z + 1+2*comput_base(i_row, j_col)
                ! if comput_base(i_row, j_col) == 0 then
                !     sum_expect_val_sigma_z = sum_expect_val_sigma_z + 1
                ! else
                !     sum_expect_val_sigma_z = sum_expect_val_sigma_z - 1
                ! end if
            end do
            long_magnetiz = long_magnetiz + evector(i_row)**2 * abs(sum_expect_val_sigma_z)
        end do
        long_magnetiz = long_magnetiz / total_sites

    end subroutine get_long_magnetiz



    subroutine get_trans_magnetiz (len_evector, evector, trans_magnetiz)
        implicit none
        integer, intent(in) :: len_evector
        real(8), intent(in) :: evector(len_evector)
        real(8), intent(out) :: trans_magnetiz

        integer :: i_row, j_col
        real(8) :: sum_trans_coeff_i

        do i_row=1, tot_Hilbert_elem
            sum_trans_coeff_i = 0
            do j_col=1, total_sites
                sum_trans_coeff_i = sum_trans_coeff_i + evector(i_row + (1+2*comput_base(i_row, j_col)) * 2**(j_col-1))
            end do
            trans_magnetiz = trans_magnetiz + evector(i_row) * sum_trans_coeff_i
        end do
        trans_magnetiz = trans_magnetiz / total_sites

    end subroutine get_trans_magnetiz    

    

    subroutine diag_lanczos(Bmat, which_eval, eval_array, evec_matrix, v0)
        implicit none

    !   %---------------------------------------%
    !   |           A*x=lambda*Bmat*x           |
    !   | where x is a vector, lambda the       |
    !   | eigenvalue to search for, and B=I.    |
    !   |                                       |
    !   | 'I' specifies that the problem to     |
    !   | solve is a standard Eigenvalues       |
    !   | problem.                              |
    !   |                                       |
    !   |              which_evals              |
    !   | 'SA' specifies that we are searching  |
    !   | for the Smallest Algebric eigenvalues.|
    !   | 'LM' for Larges Magnitude.            |
    !   %---------------------------------------%

        character,intent(in) :: Bmat(1), which_eval(2)

        real(8),intent(out) :: eval_array(nev)                            !array for eigenvalues
        real(8),intent(out),optional :: evec_matrix(tot_Hilbert_elem,nev) !1 array for each eigenvector --> matrix
        real(8),intent(in),optional :: v0(tot_Hilbert_elem)               !initial vector for lanczos
        
        integer :: stat           
        logical :: rvec = .false.   !Do you want to compute eigenvectors too?
        ido = 0                     !reverse communication variable
        tol = zero                  !if zero, machine precision is used
        stat = 0 !final state of computation. If altered an error occurred

        if (present(v0)) then 
            info=1  
            resid=v0
        else
            info=0
        end if

        if (present(evec_matrix)) rvec=.true.

        !     %-----------------------------------------------------%
        !     |                                                     |
        !     | Specification of stopping rules and initial         |
        !     | conditions before calling SSAUPD                    |
        !     |                                                     |
        !     | TOL  determines the stopping criterion.             |
        !     |                                                     |
        !     |      Expect                                         |
        !     |           abs(lambdaC - lambdaT) < TOL*abs(lambdaC) |
        !     |               computed   true                       |
        !     |                                                     |
        !     |      If TOL .le. 0,  then TOL <- macheps            |
        !     |           (machine precision) is used.              |
        !     |                                                     |
        !     | IDO  is the REVERSE COMMUNICATION parameter         |
        !     |      used to specify actions to be taken on return  |
        !     |      from SSAUPD. (See usage below.)                |
        !     |                                                     |
        !     |      It MUST initially be set to 0 before the first |
        !     |      call to SSAUPD.                                |
        !     |                                                     |
        !     | INFO on entry specifies starting vector information |
        !     |      and on return indicates error codes            |
        !     |                                                     |
        !     |      Initially, setting INFO=0 indicates that a     |
        !     |      random starting vector is requested to         |
        !     |      start the ARNOLDI iteration.  Setting INFO to  |
        !     |      a nonzero value on the initial call is used    |
        !     |      if you want to specify your own starting       |
        !     |      vector (This vector must be placed in RESID.)  |
        !     |                                                     |
        !     | The work array WORKL is used in SSAUPD as           |
        !     | workspace.  Its dimension length_workl is set as    |
        !     | illustrated below.                                  |
        !     |                                                     |
        !     %-----------------------------------------------------%
        
        !     %---------------------------------------------------%
        !     | Specification of Algorithm Mode:                  |
        !     |                                                   |
        !     | This program uses the exact shift strategy        |
        !     | (indicated by setting PARAM(1) = 1).              |
        !     | IPARAM(3) specifies the maximum number of Arnoldi |
        !     | iterations allowed.  Mode 1 of SSAUPD is used     |
        !     | (IPARAM(7) = 1). All these options can be changed |
        !     | by the user. For details see the documentation in |
        !     | SSAUPD.                                           |
        !     %---------------------------------------------------%
        
        ishfts = 1
        maxitr = 1000 !max iteration
        mode1 = 1
        
        iparam(1) = ishfts
        iparam(3) = maxitr
        iparam(7) = mode1
        
        !     %------------------------------------------------%
        !     | M A I N   L O O P (Reverse communication loop) |
        !     %------------------------------------------------%

        !        %---------------------------------------------%
        !        | Repeatedly call the routine DSAUPD and take |
        !        | actions indicated by parameter IDO until    |
        !        | either convergence is indicated or maxitr   |
        !        | has been exceeded.                          |
        !        %---------------------------------------------%

        do while (ido /= 99)

            call dsaupd ( ido, Bmat, tot_Hilbert_elem, which_eval, nev, tol, resid, ncv, &
            lancz_basis_matrix, ldv, iparam, ipntr, workd, workl, length_workl, info )

            if (ido .eq. -1 .or. ido .eq. 1) then
                
                !       %--------------------------------------%
                !       | Perform matrix vector multiplication |
                !       |              y <--- OP*x             |
                !       | The user should supply his/her own   |
                !       | matrix vector multiplication routine |
                !       | here that takes workd(ipntr(1)) as   |
                !       | the input, and return the result to  |
                !       | workd(ipntr(2)).                     |
                !       %--------------------------------------%

                stat = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.d0, sp_matrix, &
                descr_sp_matrix, workd(ipntr(1)), 0.d0, workd(ipntr(2))) 

                if (stat .ne. 0) print*, 'error mkl_sparse_d_mv: ', stat
                
                !       %-----------------------------------------%
                !       | L O O P   B A C K to call DSAUPD again. |
                !       %-----------------------------------------%
                
            end if
        end do
        
        !     %----------------------------------------%
        !     | Either we have convergence or there is |
        !     | an error.                              |
        !     %----------------------------------------%
        
        if ( info .lt. 0 ) then
            print*, 'error in diagonalization dsaupd', info
            stop
        endif
        
        !        %-------------------------------------------%
        !        | No fatal errors occurred.                 |
        !        | Post-Process using DSEUPD.                |
        !        |                                           |
        !        | Computed eigenvalues may be extracted.    |
        !        |                                           |
        !        | Eigenvectors may be also computed now if  |
        !        | desired.  (indicated by rvec = .true.)    |
        !        |                                           |
        !        | The routine DSEUPD now called to do this  |
        !        | post processing (Other modes may require  |
        !        | more complicated post processing than     |
        !        | mode1.)                                   |
        !        |                                           |
        !        %-------------------------------------------%
        
        !ALL computes NEV eigenvector, "S" computes only the ones specified by the array select
        !The firs v contains the vector B-Orthonormal in case of Generalized eigenvalues
        
        call dseupd ( rvec, 'All', select, eval_array, lancz_basis_matrix, ldv, sigma,&
        &         Bmat, tot_Hilbert_elem, which_eval, nev, tol, resid, ncv, lancz_basis_matrix, ldv,&
        &         iparam, ipntr, workd, workl, length_workl, ierr ) 
        
        !         %----------------------------------------------%
        !         | Eigenvalues are returned in the first column |
        !         | of the two dimensional array D and the       |
        !         | corresponding eigenvectors are returned in   |
        !         | the first NCONV (=IPARAM(5)) columns of the  |
        !         | two dimensional array V if requested.        |
        !         | Otherwise, an orthogonal basis for the       |
        !         | invariant subspace corresponding to the      |
        !         | eigenvalues in D is returned in V.           |
        !         %----------------------------------------------%
        
        if ( ierr .ne. 0) then
            print*, 'error in diagonalization dseupd', ierr
            stop
        endif
        
        nconv =  iparam(5)

        if (rvec) then
            do j=1,nconv
                call dcopy(tot_Hilbert_elem, lancz_basis_matrix(:,j), 1, evec_matrix(:,j), 1) !copio gli autovettori da v in evec_matrix. puoi direttamente fare assegnazione con i puntatori
            end do

        endif

    end subroutine diag_lanczos

end module Quantum_Ising_1D