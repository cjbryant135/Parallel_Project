PROGRAM GAUSS_TEST
USE MPI
USE sort_mod
use oneD_module

IMPLICIT NONE
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: x,w
INTEGER :: N, my_rank, P, ierror, i

CALL MPI_INIT(ierror)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, P, ierror)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)
IF(my_rank == 0) THEN
   WRITE(*,*) "ENTER N"
   READ(*,*) N
END IF

CALL MPI_BCAST(N, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
ALLOCATE(x(N), w(N))
CALL GaussLegendre(x,w,N,my_rank, P)

IF(my_rank == 0) THEN
   DO i = 1, N
      WRITE(*,*) w(i), x(i)
   END DO
END IF

CALL MPI_FINALIZE(ierror)

END PROGRAM GAUSS_TEST
