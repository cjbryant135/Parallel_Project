PROGRAM test_sort
USE MPI
USE sort_mod
IMPLICIT NONE
INTEGER :: my_rank, P, N, ierror, i
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: indx
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: x

CALL MPI_Init(ierror)
CALL MPI_Comm_Size(MPI_COMM_WORLD,P,ierror)
CALL MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierror)

N = 16
ALLOCATE(x(N), indx(N))
DO i = 1,N
  x(i) = COS(DBLE(i))
END DO

IF(my_rank == 0 ) THEN
  WRITE(*,*) 'before sorting:'
  DO i = 1,N
    WRITE(*,*) x(i)
  END DO
  WRITE(*,*) ''
END IF

CALL sort(x,indx,N,P,my_rank)


IF(my_rank == 0 ) THEN !something aint right here
    WRITE(*,*) 'After sorting: '
  DO i = 1,N
    WRITE(*,*) x(i)
  END DO
  WRITE(*,*)
END IF



DEALLOCATE(x,indx)

CALL MPI_Finalize(ierror)
END PROGRAM test_sort
