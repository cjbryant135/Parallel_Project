PROGRAM test
USE MPI
IMPLICIT NONE

INTEGER :: my_int, ierror, P, rank, i, j 
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: x,my_x,x_recv
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: x_two, my_x_two
INTEGER, DIMENSION(MPI_Status_Size) :: my_status

CALL MPI_Init(ierror)

CALL MPI_Comm_Size(MPI_COMM_WORLD,P,ierror)
CALL MPI_Comm_Rank(MPI_COMM_WORLD,rank,ierror)

IF(rank .EQ. 0) THEN
  ALLOCATE(x(4*P),x_two(2,4*P))
END IF

ALLOCATE(my_x(4),my_x_two(2,4),x_recv(2,1))

DO i = 1,4
  
  my_x(i) = DBLE(rank*i)
  
END DO


DO i = 1,4
  
  my_x_two(1,i) = my_x(i)
  my_x_two(2,i) = rank
END DO

my_int = 1

DO j = 1,4*P
  
  CALL MPI_Reduce(my_x_two(my_int), x_recv(1), 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, 0, &
    MPI_COMM_WORLD, ierror)
  
  WRITE(*,*) x_recv(1,1), x_recv(2,1)



END DO

DEALLOCATE(x,x_two,my_x,my_x_two, x_recv)





CALL MPI_Finalize(ierror)

END PROGRAM test
