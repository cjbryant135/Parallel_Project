PROGRAM tester
USE oneD_module
USE MPI
USE sort_mod
IMPLICIT NONE
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: x,w
INTEGER :: N=6,i

ALLOCATE(x(N),w(N))


CALL GaussLegendre(x,w,N)



!DO i = 1,4
!  WRITE(*,*) x(i)
!END DO


END PROGRAM tester

