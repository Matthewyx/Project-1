SUBROUTINE OUTPUT
   USE SIMPLE 
   IMPLICIT NONE
!!!!!!!!!!!!!!===OUTPUT===!!!!!!!!!!!!!!!!
   
OPEN(2,FILE='QRESU.PLT')
    WRITE(2,*)"VARIABLES =X, Y, C,U,DCDN"
    WRITE(2,*)"ZONE ","I=",NX+1,"J=",NY+1,",","F=BLOCK"
    DO J=0,NY
        WRITE(2,*)(X(I),I=0,NX)
    END DO
    DO J=0,NY
        WRITE(2,*)(Y(J),I=0,NX)
    END DO
    DO J=0,NY
        WRITE(2,*)(C(I,J),I=0,NX)
    END DO
    DO J=0,NY
        WRITE(2,*)(U(I,J),I=0,NX)
    END DO
    DO J=0,NY
        WRITE(2,*)(DCDN(I,J),I=0,NX)
    END DO
    CLOSE(2)
   
RETURN
END