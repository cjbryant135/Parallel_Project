MODULE oneD_module
USE MPI
USE sort_mod
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

IF(my_rank .NE. 0) THEN
   DEALLOCATE(x, indx)
ELSE
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



END MODULE oneD_module

