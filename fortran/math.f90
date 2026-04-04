module math
    use parkind,   only:im=>kind_im, rb=>kind_rb
    implicit none
    public :: inv
    contains

    function inv(A) result(Ainv)
        ! Returns the inverse of a matrix calculated by finding the LU
        ! decomposition.  Depends on LAPACK.

        real(kind=rb), dimension(:,:), intent(in) :: A
        real(kind=rb), dimension(size(A,1),size(A,2)) :: Ainv

        real(kind=rb), dimension(size(A,1)) :: work  ! work array for LAPACK
        integer(kind=im), dimension(size(A,1)) :: ipiv   ! pivot indices
        integer(kind=im) :: n, info

        ! External procedures defined in LAPACK
        external DGETRF
        external DGETRI

        ! Store A in Ainv to prevent it from being overwritten by
        ! LAPACK
        Ainv = A
        n = size(A,1)

        ! DGETRF computes an LU factorization of a general
        ! M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call DGETRF(n, n, Ainv, n, ipiv, info)

        if (info /= 0) then
            stop 'Matrix is numerically singular!'
        end if

        ! DGETRI computes the inverse of
        ! a matrix using the LU
        ! factorization
        ! computed by DGETRF.
        call DGETRI(n, Ainv, n, ipiv, work, n, info)

        if (info /= 0) then
            stop 'Matrix inversion failed!'
        end if
    end function inv

end module
