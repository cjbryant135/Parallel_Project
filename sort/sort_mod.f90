MODULE sort_mod
USE MPI
SUBROUTINE sort(x,indx,N,P,my_rank)
IMPLICIT NONE
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: x,indx,tempx,tempindx
INTEGER :: N,P,my_rank,leftover,div,i,j,tempx_size,temp_rank,merge_count

leftover = MOD(N,P)
div = N - leftover

IF(my_rank .GE. leftover) THEN
  ALLOCATE(tempx(div+1),tempindx(div+1))
  tempx_size = div+1
ELSE
  ALLOCATE(tempx(div),tempindx(div))
  tempx_size = div
END IF

DO i = 1+my_rank*div, my_rank*div+div
  DO j = 1, div
    tempx(j) = x(i)
    tempindx(j) = i
  END DO
END DO

IF(my_rank .GE. leftover) THEN
  tempx(div+1) = x(P*div + my_rank-leftover+1)
  tempindx(div+1) = P*div + my_rank-leftover+1
END IF

CALL MergeSort(tempx,tempx_size,tempindx)

!Put it all back together!





END DO








DEALLOCATE(tempx,tempindx)
END SUBROUTINE sort


RECURSIVE SUBROUTINE MergeSort(tempx,tempx_size,tempindx)
IMPLICIT NONE
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: tempx, tempindx, left, right, leftindx, rightindx
INTEGER tempx_size,left_size,right_size, pivot, i

IF(tempx_size .LE. 1) THEN 
  RETURN
END IF

pivot = FLOOR(tempx_size / 2)

IF(MOD(tempx_size,2) .EQ. 0) THEN
  ALLOCATE(left(pivot),right(pivot),leftindx(pivot),rightindx(pivot))
  left_size = pivot
  right_size = pivot
ELSE
  ALLOCATE(left(pivot),right(pivot+1),leftindx(pivot),rightindx(pivot+1))
  left_size = pivot
  right_size = pivot+1
END IF

DO i = 1,pivot
  left(i) = tempx(i)
  leftindx(i) = tempindx(i)
  right(i) = tempx(pivot+i)
  rightindx(i) = tempindx(pivot+i)
END DO

IF(MOD(tempx_size,2) .NE. 0) THEN 
  right(tempx_size) = tempx(tempx_size)
  rightindx(tempx_size) = tempindx(tempx_size)
END IF

CALL MergeSort(left,left_size,leftindx)
CALL MergeSort(right,right_size,rightindx)

CALL MergeIt(left,left_size,right,right_size,leftindx,rightindx,tempx,tempindx)

DEALLOCATE(left,right,leftindx,rightindx)
END RECURSIVE SUBROUTINE MergeSort

SUBROUTINE MergeIt(left,left_size,right,right_size,leftindx,rightindx,tempx,tempindx)
IMPLICIT NONE
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: tempx, tempindx, left, right, leftindx, rightindx
INTEGER tempx_size,left_size,right_size, i, templ, tempr, k

templ = 1
tempr = 1
k = 1

WHILE(templ .LE. left_size) DO
  WHILE(tempr .LE. right_size) DO
    
    IF(left(templ) .GE. right(tempr)) THEN
      
      tempx(k) = right(tempr)
      tempindx(k) = rightindx(tempr)
      k = k+1
      tempr = tempr+1
    ELSE
      tempx(k) = left(templ)
      tempindx(k) = leftindx(templ)
      k = k+1
      templ = templ+1

    END IF
    
    IF(tempr .GT. right_size) THEN
      
      DO i = templ, left_size
        tempx(k) = left(i)
        tempindx(k) = leftindx(i)
        k = k+1
      END DO
      
      templ = left_size+1

    END IF
    
    
    IF(templ .GT. left_size) THEN
      
      DO i = tempr, right_size
        tempx(k) = right(i)
        tempindx(k) = rightindx(i)
        k = k+1
      END DO
      
      tempr = right_size+1

    END IF
    
  END DO


END DO

END SUBROUTINE MergeIt





END MODULE sort_mod
