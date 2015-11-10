MODULE oneD_module
USE MPI
USE sort_mod
USE oned_parameters
CONTAINS

SUBROUTINE GaussLegendre(x, w, N, my_rank, P)
  IMPLICIT NONE
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: D,x,w,beta,Lambda,Work,indx
  INTEGER :: N,i,j,Lwork,Info,num_eig,Liwork, my_rank, P, ierror
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: V, T, Eye
  DOUBLE PRECISION :: abs_tol
  INTEGER, ALLOCATABLE, DIMENSION(:) :: Isuppz,Iwork
 
 !BEGIN PROGRAM HERE
  ALLOCATE(indx(N))
  IF(my_rank == 0) THEN
  
     ALLOCATE(beta(N),T(N,N))
  
  T(1:N,1:N) = 0
  !TEST
  !  DO i = 1, N
  !   DO j=1,N
  !    WRITE(*,*) T(i,j)
  !  END DO
 ! END DO
  DO i=1,N-1
    beta(i) = 0.5d0 / SQRT(1.d0 - ( 2.d0 * DBLE(i))**(-2.d0))
    T(i,i+1) = beta(i)
    T(i+1,i) = beta(i)    
  END DO
  !TEST
  !  DO i = 1, N
  !   DO j=1,N
  !     WRITE(*,*) T(i,j)
  !   END DO
  ! END DO

  Lwork = 18*N+1
  Liwork = 10*N+1
  ALLOCATE(V(N,N), D(N), Lambda(N),Isuppz(2*N), Work(Lwork),&
    Iwork(Liwork))
  !MAKE IDENTITY MATRIX
  
  ALLOCATE(Eye(N,N))
  Eye(1:N,1:N) = 0
  DO i = 1,N
    Eye(i,i) = 1
    D(i) = T(i,i)
  END DO
  
  
  
  !CALL SGGEV('N','V',N,T,N,Eye,N,Lambda_num_real,Lambda_num_im,Lambda_denom_real,V,N,V,N,Work,Lwork,Info)  
  !DOES NOT WORK YET!
  
  

  CALL DSTEGR('V','A',N, D, beta,0.d0,0.d0,0,0,abs_tol,num_eig,Lambda,V,N,Isuppz,Work,Lwork,Iwork,&
    Liwork,Info) !WORKING!

  !
  !
  DO i = 1, N
    x(i) = Lambda(i) !put eigenvalues into x
  END DO

END IF
  
CALL MPI_Bcast(x, N, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
DO i = 1, N
   indx(i) = DBLE(i)
END DO

CALL sort(x, indx, N, P, my_rank)

IF(my_rank == 0) THEN
   DO i = 1, N
      Lambda(i) = x(i)
      w(i) = DBLE(2)*V(1,INT(indx(i)))**2
   END DO
END IF


  !TEST STUFF
  !DO i = 1, N
  !  DO j=1,N
  !    WRITE(*,*) T(i,j)
  !  END DO
  !END DO
  !Stores eigenvalues in x
  !DO i = 1,N
  !  WRITE(*,*)  V(1:N,i)
    
  !END DO
  
  !NEED TO SORT EIGENVALUES AND TRACK CHANGES IN index
  


IF(my_rank == 0) THEN
   DEALLOCATE(beta,T,V,D,Lambda,Isuppz,Work,Iwork)
END IF
   !DEALLOCATE(beta) 
   !WRITE(*,*) 'beta'
   !DEALLOCATE(T)
   !WRITE(*,*) 'T'
   !DEALLOCATE(V)
   !WRITE(*,*) 'V'
   !DEALLOCATE(D)
   !WRITE(*,*) 'D'
   !DEALLOCATE(Lambda)
   !WRITE(*,*) 'Lambda'
   !DEALLOCATE(Isuppz)
   !WRITE(*,*) 'Isuppz'
   !DEALLOCATE(Work)
   !WRITE(*,*) 'Work'
   !DEALLOCATE(Iwork)
   !WRITE(*,*) 'Iwork'



END SUBROUTINE GaussLegendre


SUBROUTINE RTE_oneD(N, my_rank, P, Intensity) !We assume intensity is already allocated, size N
!x w N my_rank P
IMPLICIT NONE
INTEGER :: P, my_rank, N, i, j, LWORK, Info, ierror
INTEGER, PARAMETER :: master = 0
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: x, w, lambda_num, lambda_de, alphai, Work, &
  lambda, evals, indx, cond, ab, Intensity
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: B, HH, A, eye, VL, VR, V, cutV, U, Uplus, Uminus, Vplus, Vminus, block, Psi
INTEGER, ALLOCATABLE, DIMENSION(:) :: IPIV


ALLOCATE(x(N), w(N), lambda(N), indx(N))
CALL GaussLegendre(x, w, N, my_rank, P)

IF(my_rank == master) THEN
  ALLOCATE(B(N,N), HH(N,N), A(N,N), eye(N,N))
  B(:,:) = 0
  LWORK = MAX(1,8*N)+1
  ALLOCATE(Work(LWORK))
  ALLOCATE(lambda_num(N), lambda_de(N), alphai(N))
  ALLOCATE(VL(N,N), VR(N,N))
  DO i = 1, N
    B(i,i) = x(i)
    eye(i,i) = 1.d0
    DO j = 1,N
      HH(i,j) = h(x(i),x(j))*w(j)
    END DO
  END DO
  
  !DO i = 1, N
  !   DO j = 1, N
  !      WRITE(*,*) B(i,j)
  !   END DO
  !END DO
  
  A(:,:) = -(eye(:,:)-albedo()*HH(:,:))
 

  CALL DGEGV('N', 'V', N, A, N, B, N, lambda_num, alphai, lambda_de, VL, N, VR, N, Work, LWORK, Info) 
  lambda(:) = lambda_num(:) / lambda_de(:)
  !WRITE(*,*) lambda
  DEALLOCATE(Work, VL, alphai, A, B, HH, eye, lambda_num, lambda_de)
  
ELSE
  DEALLOCATE(x,w)
END IF

CALL MPI_Bcast(lambda, N, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierror) !everyone has lambda

DO i = 1, N
  indx(i) = DBLE(i) !create indx vector
END DO


CALL sort(lambda, indx, N, P, my_rank) !eigenvalues sorted and changes tracked in indx
IF(my_rank == master) THEN
  !WRITE(*,*) indx
  ALLOCATE(evals(N/2))
  evals(:) = lambda(N/2+1:N) !store only positive eigenvalues
  DEALLOCATE(lambda)
  ALLOCATE(V(N,N))
  DO i = 1, N
    V(1:N,i) = VR(1:N,INT(indx(i)))
  END DO

  !WRITE(*,*) "V"
  !DO i = 1, N
  !   DO j = 1, N
  !      WRITE(*,*) V(i,j)
  !   END DO
  !END DO

  DEALLOCATE(VR)
  ALLOCATE(cutV(N,N/2), U(N,N/2))
  cutV = V(:,N/2+1:N)
  DEALLOCATE(V)

  CALL flipud(cutV, U, N)
  !DO i = 1, N
  !   DO J = 1, N/2
  !      WRITE(*,*) U(i,j)
  !   END DO
  !END DO
  !WRITE(*,*) ""

  ALLOCATE(Uplus(N/2,N/2), Uminus(N/2,N/2), Vplus(N/2,N/2), Vminus(N/2,N/2))
  !Construct blocks to put into block matrix
  Uminus = U(1:N/2,:)
  Uplus = U(N/2+1:N,:)

  Vplus = cutV(N/2+1:N,:)
  Vminus = cutV(1:N/2,:)
  

  ALLOCATE(block(N,N))
  !set the diagonal blocks of block
  block(1:N/2, 1:N/2) = Uplus(:,:)
  block(N/2+1:N, N/2+1:N) = Vminus(:,:)
  
  !set psi up
  ALLOCATE(Psi(N/2, N/2)) 
  Psi(:,:) = 0
  DO i = 1, N/2
    Psi(i,i) = EXP(-evals(i)*taufinal())
  END DO
  
  block(1:N/2,N/2+1:N) = MATMUL(Vplus, Psi) 
  block(N/2+1:N, 1:N/2) = MATMUL(Uminus,Psi)

  !DO i = 1, N
  !   DO j = 1, N
  !      WRITE(*,*) block(i,j)
  !   END DO
  !END DO
  
  DEALLOCATE(Uplus, Uminus, Vplus, Vminus, Psi)
  ALLOCATE(cond(N))
  cond(1:N/2) = alpha(x(N/2+1:N), N/2)
  cond(N/2+1:N) = beta(x(1:N/2), N/2)
  

  
  !NOW SOLVE THAT SYSTEM 
  ALLOCATE(IPIV(N))
  CALL DGESV(N, 1, block, N, IPIV, cond, N, Info)
  ALLOCATE(ab(N))
  ab = cond

  !DO i = 1, N
  !   WRITE(*,*) cond(i)
  !END DO
  WRITE(*,*) ""
 
  DEALLOCATE(block, cond, x, w, IPIV)
  
  Intensity(:) = 0.d0

  DO i = 1, N/2
    Intensity(:) = Intensity(:) + U(:,i)*EXP(-evals(i)*tau())*ab(i) + &
      cutV(:,i)*EXP(evals(i)*(tau() - taufinal() ))*ab(N/2+i)
  END DO
  


ELSE 
  DEALLOCATE(lambda, indx)
END IF




END SUBROUTINE RTE_oneD


SUBROUTINE flipud(cutV, U, N)
IMPLICIT NONE
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: cutV, U
INTEGER :: i, N

DO i = 1, N
  U(i,:) = cutV(N-i+1,:)
END DO

END SUBROUTINE










END MODULE oneD_module

