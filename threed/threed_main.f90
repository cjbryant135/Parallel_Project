PROGRAM threed_main
USE MPI
USE threed_mod
USE threed_parameters
IMPLICIT NONE

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: k, mu, phi, w
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: my_PP, my_eye, my_Mu, my_Phi
INTEGER :: P, P2, my_rank, my_rank2, i, j, l, left_over, new_comm, color_id, div, k1_size, &
  ierror, fail, low, high, my_r, left_PP, div_PP, row_PP
INTEGER, PARAMETER :: master = 0
DOUBLE PRECISION :: dk, eta, k1, k2
DOUBLE PRECISION, DIMENSION(2) :: times

CALL MPI_INIT(ierror)
WRITE(*,*) 'init'
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, P, ierror)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)
times(1) = MPI_Wtime()
WRITE(*,*) 'I started'

!Parameters
!x_max = 1d0
!N = 100
!pi = 3.141592653589793d0
dk = pi()/x_max()
!master = 0
!N_mu = 6
!N_phi = 4


!P MUST BE PERFECT SQUARE, CHECK THIS CONDITION, NOTE ALSO N > SQRT(P)!
eta = SQRT(DBLE(P))
!IF(ABS(eta - FLOOR(eta)) > 0.0001 .OR. N() < eta) THEN
!  WRITE(*,*) 'Problem with N or eta size!!!!'
!  CALL MPI_Abort(MPI_COMM_WORLD, fail, ierror)
!END IF

WRITE(*,*) 'I didnt abort'

!BECAUSE WE ASSUME A SQUARE GRID!
ALLOCATE( k(N()) )

IF (my_rank == master) THEN
   DO i = 1, N()
      k(i) = -N()/2 + i-1 
   END DO

   k = k*dk
END IF


ALLOCATE(mu(N_mu()), w(N_mu()), phi(N_phi()))
CALL GaussLegendre(mu, w, N_mu(), my_rank, P) !get mu and w
  
IF(my_rank == master) THEN
  !create phi
  DO i = 0, N_phi() - 1
    phi(i+1) = (DBLE(2)*pi()/DBLE(N_phi()))*DBLE(i)
  END DO
  
END IF

WRITE(*,*) 'I got gauss'

CALL MPI_BCAST(mu, N_mu(), MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierror)

CALL MPI_BCAST(w, N_mu(), MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierror)

CALL MPI_BCAST(phi, N_phi(), MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierror)

color_id = MOD(my_rank, INT(eta))
CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, color_id, my_rank, new_comm, ierror)
CALL MPI_COMM_SIZE(new_comm, P2, ierror)
CALL MPI_COMM_RANK(new_comm, my_rank2, ierror)

left_PP = MOD(N_mu()*N_phi(), P2)
div_PP = (N_mu()*N_phi()-left_PP)/P2

IF(my_rank2 .GE. P2-left_PP) THEN !IF you're getting extra
  my_r = div_PP + 1
ELSE
  my_r = div_PP
END IF

ALLOCATE(my_PP(my_r*N_phi(), N_mu()*N_phi()), my_eye(my_r*N_phi(), N_mu()*N_phi()), &
  my_Phi(my_r*N_phi(), N_mu()*N_phi()), my_Mu(my_r*N_phi(), N_mu()*N_phi()))

!TEST CODE
!WRITE(*,*) 'My rank is:', my_rank, ' my_r is ', my_r
!WRITE(*,*) ''
!GETTING CORRECT my_r values
!IF(my_rank == master) THEN
!  WRITE(*,*) 'div_PP = ', div_PP 
!  WRITE(*,*) ''
!END IF
!GETTING CORRECT div_P(:,:) = 0.d0
my_eye(:,:) = 0.d0
my_Phi(:,:) = 0.d0
my_Mu(:,:) = 0.d0
my_PP(:,:) = 0.d0

DO i = 1, my_r !Create my matrices
  !MU CONSTRUCTION IS WORKING 
  IF(my_r == div_PP) THEN
    my_Mu(i,my_rank2*div_PP+i) = mu(INT(CEILING(DBLE(my_rank2*div_PP+i)/DBLE(N_phi()))))
  ELSE 
    my_Mu(i,my_rank2*div_PP+i+(my_rank2-(P2-left_PP))) = & 
      mu(INT(CEILING(DBLE(my_rank2*div_PP+i+(my_rank2-(P2-left_PP)))/DBLE(N_phi()))))
  END IF
  
  !PHI CONSTRUCTION
  IF(my_r == div_PP) THEN
    my_Phi(i,my_rank2*div_PP+i) = phi(MOD(my_rank2*div_PP+i-1, N_phi())+1)
  ELSE
    my_Phi(i,my_rank2*div_PP+i+(my_rank2-(P2-left_PP))) = &
      phi(MOD(my_rank2*div_PP+i+(my_rank2-(P2-left_PP))-1, N_phi())+1)
  END IF

  !IDENTITY CONSTRUCTION
  IF(my_r == div_PP) THEN
    my_eye(i, my_rank2*div_PP+i) = 1.d0
  ELSE 
    my_eye(i, my_rank2*div_PP+i+(my_rank2-(P2-left_PP))) = 1.d0
  END IF

  !P construction :( 
  IF(my_r == div_PP) THEN
    row_PP = my_rank2*div_PP+i
  ELSE
    row_PP = my_rank2*div_PP+i+(my_rank2-(P2-left_PP))
  END IF
  DO j = 1, N_mu()*N_phi() !do across this entire row
    CALL PP_entry(row_PP, j, mu, phi, w, my_PP, i)
  END DO

END DO

!IF(my_rank == 1) THEN
!  WRITE(*,*) 'for core 1, my_Mu(1,7) = ', my_Mu(1,7)
!END IF
!WRITE(*,*) 'eta is ', eta
!TEST CODE FOR P
!IF(color_id == 0) THEN
!  DO i = 0, INT(eta) -1
!    CALL MPI_Barrier(new_comm, ierror)
!    IF(my_rank2 == i) THEN
!      WRITE(*,*) 'I am core: ', my_rank2
!      WRITE(*,*) '-------------------------------------------------------------'
!      DO j = 1, my_r
!        WRITE(*,*) my_PP(j,:)
!        WRITE(*,*) ''
!      END DO
!      WRITE(*,*) '--------------------------------------------------------------'
!    END IF
!
!
!  END DO
!
!
!END IF

times(2) = MPI_Wtime()

IF(my_rank == master) THEN
  WRITE(*,*) 'TIME! = ', times(2) - times(1)
END IF

DEALLOCATE(mu, w, phi, my_eye, my_PP, my_Mu, my_Phi)












!Test code for my_Mui
!DO i = 0, P-1 
!  CALL MPI_Barrier(new_comm, ierror)
!  IF(my_rank2 == master .AND. i == 0) THEN
!    WRITE(*,*) ''
!    WRITE(*,*) 'Mu is: '
!    WRITE(*,*) 'Phi is: '
!    DO j = 1, N_phi()
!      WRITE(*,*) mu(j)
!       WRITE(*,*) phi(j)
!    END DO
!    WRITE(*,*) ''
!  END IF
!  IF(my_rank2 == i) THEN
!    WRITE(*,*) 'I am core: ', my_rank2
!    DO j = 1, my_r
!      IF(my_r == div_PP) THEN
!        WRITE(*,*) my_Phi(j, my_rank2*div_PP+j)
!      ELSE
!        WRITE(*,*) my_Phi(j, my_rank2*div_PP+j+(my_rank2-(P2-left_PP)))
!      END IF
!    END DO
!    WRITE(*,*) ''
!  END IF
!END DO

!DO i = 0, P-1 
!  CALL MPI_Barrier(new_comm, ierror)
!  IF(my_rank2 == i) THEN
!    WRITE(*,*) 'I am: ', my_rank2
!    WRITE(*,*) '------------------------------------------------------------------------'
!    DO j = 1, my_r
!      WRITE(*,*) my_Mu(j,:)
!      WRITE(*,*) ''
!    END DO
!    WRITE(*,*) '------------------------------------------------------------------------'
!  END IF
!END DO
















!SHIT WITH K below here
CALL MPI_BCAST(k, N(), MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierror)



!Communicator is split yo

left_over = MOD(N(), INT(eta))
div = (N() - left_over) / eta

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
DO i = 1, N()
  k2 = k(i) 
  DO j = low, high
    
    CALL MPI_BARRIER(new_comm, ierror)

    k1 = k(j)
    










































  END DO
END DO





DEALLOCATE(k)






CALL MPI_FINALIZE(ierror)




END PROGRAM
