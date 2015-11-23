PROGRAM scala_test
USE MPI 
IMPLICIT NONE
 
INTEGER :: P, my_rank, ierror, i, j, INFO
INTEGER, PARAMETER :: master = 0, LWORK = 1000
INTEGER, DIMENSION(4) :: DESCA
DOUBLE PRECISION, DIMENSION(4) :: my_A
DOUBLE PRECISION, DIMENSION(2) :: TAU
DOUBLE PRECISION, DIMENSION(LWORK) :: WORK
CALL MPI_Init(ierror)
CALL MPI_COMM_Size(MPI_COMM_WORLD, P, ierror)
CALL MPI_COMM_Rank(MPI_COMM_WORLD, my_rank, ierror)

CALL DESCINIT(DESCA, 4, 4, 1, 4, 0, 0, my_rank, 1, INFO)
!What it means
!DESCINIT(DESCA, number of rows of A, cols of A, rows of my_A, cols of my_A, 0, 0, my_rank, rows of my_A, info) 


!Create A
IF(my_rank == master) THEN
   my_A = [1, 4, 6, -2]
ELSE IF(my_rank == 1) THEN
   my_A = [-2, 5, 1, 0]
ELSE IF(my_rank == 2) THEN
  my_A = [0, 3, 6, 3] 
ELSE
  my_A = [2, 3, 2, -1]
END IF

CALL PDGEHRD(4, 1, 4, my_A, my_rank+1, 1, DESCA, TAU, WORK, LWORK, INFO)
!PDGEHRD(N, ?, ?, my_A, first row of my_A in A, first col of my_A in A, DESCA, 





CALL MPI_Finalize(ierror)




END PROGRAM scala_test
