PROGRAM threed_main
USE MPI
USE threed_mod
IMPLICIT NONE

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: k
INTEGER :: N, P, P2, my_rank, my_rank2, i, left_over, new_comm, color_id, div, k1_size, master, ierror
DOUBLE PRECISION :: pi, x_max, y_max, dk, eta

CALL MPI_INIT(ierror)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, P, ierror)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)

x_max = 1d0
pi = 3.141592653589793d0
dk = pi/x_max
master = 0
!P MUST BE PERFECT SQUARE, CHECK THIS CONDITION, NOTE ALSO N > SQRT(P)!
eta = sqrt(DBLE(P))

!BECAUSE WE ASSUME A SQUARE GRID!
ALLOCATE( k(N) )

IF (my_rank == master) THEN
   DO i = 1, N
      k(i) = -N/2 + i-1 
   END DO

   k = k*dk
END IF

CALL MPI_BCAST(k, N, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierror)

color_id = mod(my_rank, int(eta))
CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, color_id, my_rank, new_comm, ierror)

CALL MPI_FINALIZE(ierror)




END PROGRAM
