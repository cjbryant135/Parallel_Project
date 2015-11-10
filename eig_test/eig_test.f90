PROGRAM eig_test

IMPLICIT NONE
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: A, B, VL, VR
INTEGER :: N, i, j, LWORK, Info
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: alphar, alphai, beta, Work

WRITE(*,*) 'Enter a value for N: '
READ(*,*) N

LWORK = MAX(1,8*N) +1

ALLOCATE(Work(LWORK))
ALLOCATE(A(N,N), B(N,N), alphar(N), alphai(N), beta(N))

ALLOCATE(VL(N,N))
ALLOCATE(VR(N,N))

B(:,:) = 0

DO i = 1, N
  DO j = 1, N
    A(i,j) = DBLE(i)**2+3d0*DBLE(j)+DBLE(i)**j
  END DO
  !A(i,i) = 1
  IF(i .LE. N/2) THEN
    B(i,i) = -1
  ELSE
    B(i,i) = 1
  END IF
END DO




!solve for eigs
CALL DGEGV('N', 'V', N, A, N, B, N, alphar, alphai, beta, VL, N, VR, N, Work, LWORK, Info)   

DO i = 1,N
  WRITE(*,*) VR(i,1), VR(i,2), VR(i,3), VR(i,4)
  !WRITE(*,*) alphar(i)/beta(i)
END DO



DEALLOCATE(A, B, alphar, alphai, beta, VL, VR, Work)

END PROGRAM eig_test
