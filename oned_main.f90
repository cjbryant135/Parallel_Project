PROGRAM oned_main
USE MPI
USE oned_module

IMPLICIT NONE
INTEGER :: ierror, P, my_rank, N, i
DOUBLE PRECISION :: endTime, startTime
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: Intensity, x, w, lambda, indx
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: VR
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

IF(my_rank == master) THEN
  ALLOCATE(Intensity(N))
END IF

ALLOCATE(x(N), w(N), lambda(N), indx(N))

CALL GaussLegendre(x, w, N, my_rank, P)


IF(my_rank == master) THEN
  ALLOCATE(VR(N,N))
  CALL eig_solve(x, w, N, lambda, VR)

END IF

CALL MPI_Bcast(lambda, N, MPI_DOUBLE_PRECISION, master, &
  MPI_COMM_WORLD, ierror)

DO i = 1, N
  indx(i) = DBLE(i)
END DO

!IF(my_rank == master) THEN
!  WRITE(*,*) 'Before sorting '
!  DO i = 1, N
!    WRITE(*,*) lambda(i), indx(i)
!  END DO
!  WRITE(*,*) '---------------------'
!END IF

CALL sort(lambda, indx, N, P, my_rank)

!IF(my_rank == master) THEN
!  WRITE(*,*) 'After sorting'
!  DO i = 1, N
!    WRITE(*,*) lambda(i), indx(i)
!  END DO
!  WRITE(*,*) ''
!END IF




IF(my_rank == master) THEN
  
  CALL get_Intensity(lambda, VR, indx, x, w, N, Intensity)
  WRITE(*,*) ''
  !WRITE(*,*) 'Intensity' 
  !DO i = 1, N
  !  WRITE(*,*) Intensity(i)
  !END DO
  

  DEALLOCATE(Intensity, VR)
END IF



DEALLOCATE(lambda, indx, x, w)



endTime = MPI_wtime()
IF(my_rank == 0) THEN
WRITE(*,*) endTime - startTime, "seconds"
END IF

CALL MPI_Finalize(ierror)

END PROGRAm oned_main
