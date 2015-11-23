PROGRAM threed_main
USE MPI
USE threed_mod
IMPLICIT NONE

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: k, mu, phi, w
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: my_PP, my_eye, my_Mu, my_Phi
INTEGER :: N, P, P2, my_rank, my_rank2, i, j, left_over, new_comm, color_id, div, k1_size, &
  master, ierror, fail, low, high, N_mu, N_phi, my_r, left_PP, div_PP
DOUBLE PRECISION :: pi, x_max, y_max, dk, eta, k1, k2

CALL MPI_INIT(ierror)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, P, ierror)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)

!Parameters
x_max = 1d0
N = 100
pi = 3.141592653589793d0
dk = pi/x_max
master = 0
N_mu = 16
N_phi = 8


!P MUST BE PERFECT SQUARE, CHECK THIS CONDITION, NOTE ALSO N > SQRT(P)!
eta = SQRT(DBLE(P))
IF(ABS(eta - FLOOR(eta)) > 0.0001 .OR. N < eta) THEN
  WRITE(*,*) 'Problem with N or eta size!!!!'
  CALL MPI_Abort(MPI_COMM_WORLD, fail, ierror)
END IF



!BECAUSE WE ASSUME A SQUARE GRID!
ALLOCATE( k(N) )

IF (my_rank == master) THEN
   DO i = 1, N
      k(i) = -N/2 + i-1 
   END DO

   k = k*dk
END IF


ALLOCATE(mu(N_mu), w(N_mu), phi(N_phi))
CALL GaussLegendre(mu, w, N_mu, my_rank, P) !get mu and w
  
IF(my_rank == master) THEN
  !create phi
  DO i = 0, N_phi - 1
    phi(i+1) = (DBLE(2)*pi/DBLE(N_phi))*DBLE(i)
  END DO
  
END IF



CALL MPI_BCAST(mu, N_mu, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierror)

CALL MPI_BCAST(w, N_mu, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierror)

CALL MPI_BCAST(phi, N_phi, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierror)

left_PP = MOD(N_mu*N_phi, P)
div_PP = (N_mu*N_phi-left_PP)/P

IF(my_rank .GE. P-left_PP) THEN !IF you're getting extra
  my_r = div_PP + 1
ELSE
  my_r = div_PP
END IF

ALLOCATE(my_PP(my_r*N_phi, N_mu*N_phi), my_eye(my_r*N_phi, N_mu*N_phi), &
  my_Phi(my_r*N_phi, N_mu*N_phi), my_Mu(my_r*N_phi, N_mu*N_phi))

my_PP(:,:) = 0.d0
my_eye(:,:) = 0.d0
my_Phi(:,:) = 0.d0
my_Mu(:,:) = 0.d0











!SHIT WITH K below here
CALL MPI_BCAST(k, N, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierror)

color_id = MOD(my_rank, INT(eta))
CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, color_id, my_rank, new_comm, ierror)
CALL MPI_COMM_SIZE(new_comm, P2, ierror)
CALL MPI_COMM_RANK(new_comm, my_rank2, ierror)


!Communicator is split yo

left_over = MOD(N, INT(eta))
div = (N - left_over) / eta

IF(my_rank2 == 0) THEN !All submasters
  IF(color_id .GE. P-left_over) THEN
    k1_size = div + 1
    low = color_id*div+(color_id-1)
    high = (color_id+1)*div+(color_id-1)
  ELSE
    k1_size = div
    low = color_id*div+1
    high = (color_id+1)*div
  END IF
   
END IF

CALL MPI_BCAST(k1_size, 1, MPI_INTEGER, master, new_comm, ierror)
CALL MPI_BCAST(low, 1, MPI_INTEGER, master, new_comm, ierror)
CALL MPI_BCAST(high, 1, MPI_INTEGER, master, new_comm, ierror)








!each communicator to loop over a subset of the k combos
DO i = 1, N
  k2 = k(i) 
  DO j = low, high
    
    CALL MPI_BARRIER(new_comm, ierror)

    k1 = k(j)
    










































  END DO
END DO








DEALLOCATE(mu, w, phi, k)







CALL MPI_FINALIZE(ierror)




END PROGRAM
