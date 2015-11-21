PROGRAM one

USE MPI
IMPLICIT NONE

INTEGER ::  ierror, numCores, numCores2, myRank,  myRank2, newComm, colorId, rankId

CALL MPI_INIT(ierror)
CALL MPI_COMM_SIZE (MPI_COMM_WORLD, numCores, ierror)
CALL MPI_COMM_RANK (MPI_COMM_WORLD, myRank, ierror)

WRITE(*,*) "I AM PROCESSOR NUMBER" , myRank, "IN MPI_COMM_WORLD OF SIZE", numcores

colorId =  MOD(myRank+1, 2)
rankId = myRank

CALL MPI_COMM_split(MPI_COMM_WORLD, colorId, rankId, newComm, ierror)

CALL MPI_COMM_SIZE(newComm, numCores2, ierror)
CALL MPI_COMM_RANK(newComm, myRank2, ierror)
WRITE(*,*) "I AM PROCESSOR NUMBER", myRank2, "In newComm of size", numCores2

CALL MPI_FINALIZE(ierror)



END PROGRAM one
