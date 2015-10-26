MODULE oneD_module
USE MPI
CONTAINS

SUBROUTINE GaussLegendre(x, w, N)
  IMPLICIT NONE
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: x,w,beta,Lambda_num_real, Lambda_denom_real,Lambda_num_im, Work
  INTEGER :: N,i,j,Lwork,Info
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: V, D, T, Eye
  
  !BEGIN PROGRAM HERE
  ALLOCATE(beta(N))
  ALLOCATE(T(N,N))
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

  !MAKE IDENTITY MATRIX
  ALLOCATE(Eye(N,N))
  Eye(1:N,1:N) = 0
  DO i = 1,N
    Eye(i,i) = 1
  END DO
  
  Lwork = -1

  ALLOCATE(V(N,N), D(N,N), Lambda_num_real(N), Lambda_denom_real(N), Lambda_num_im(N), Work(1))
  CALL SGEGV('N','V',N,T,N,Eye,N,Lambda_num_real,Lambda_num_im,Lambda_denom_real,V,N,V,N,Work,Lwork,Info)  
  !DOES NOT WORK YET!

  !Stores eigenvalues in x
  DO i = 1,N
    x(i) = Lambda_num_real(i)/Lambda_denom_real(i)
  END DO
  
  !NEED TO SORT EIGENVALUES AND TRACK CHANGES IN index




  DEALLOCATE(beta,T)


END SUBROUTINE GaussLegendre


SUBROUTINE sort(x,indx,N)
IMPLICIT NONE
INTEGER, ALLOCATABLE, DIMENSION(:) :: indx  
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: x
INTEGER :: N,i








END SUBROUTINE sort
























END MODULE oneD_module

