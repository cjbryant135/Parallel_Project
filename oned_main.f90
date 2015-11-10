PROGRAM oned_main
USE MPI
USE oned_module

IMPLICIT NONE
INTEGER :: ierror, P, my_rank, N, i
DOUBLE PRECISION :: endTime, startTime
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: Intensity
INTEGER, PARAMETER :: master = 0

CALL MPI_Init(ierror)
CALL MPI_Comm_Size(MPI_COMM_WORLD, P, ierror)
CALL MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierror)

IF(my_rank == master) THEN
  WRITE(*,*) 'ENTER A VALUE FOR N: '
  READ(*,*) N
END IF
startTime = MPI_wtime()
CALL MPI_Bcast(N, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierror)

ALLOCATE(Intensity(N))

CALL RTE_oneD(N, my_rank, P, Intensity)

!test code
!IF(my_rank == master) THEN
  !DO i = 1, N
    !WRITE(*,*) Intensity(i)
  !END DO
!END IF

DEALLOCATE(Intensity)
endTime = MPI_wtime()
IF(my_rank == 0) THEN
WRITE(*,*) endTime - startTime, "seconds"
END IF

CALL MPI_Finalize(ierror)

END PROGRAm oned_main
