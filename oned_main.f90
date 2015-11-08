PROGRAM oned_main
USE MPI
USE oned_module

IMPLICIT NONE
INTEGER :: ierror, P, my_rank, N
INTEGER, PARAMETER :: master = 0

CALL MPI_Init(ierror)
CALL MPI_Comm_Size(MPI_COMM_WORLD, P, ierror)
CALL MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierror)

IF(my_rank == master) THEN
  WRITE(*,*) 'ENTER A VALUE FOR N: '
  READ(*,*) N
END IF

CALL MPI_Bcast(N, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierror)

CALL RTE_oneD(N, my_rank, P)



CALL MPI_Finalize(ierror)

END PROGRAm oned_main
