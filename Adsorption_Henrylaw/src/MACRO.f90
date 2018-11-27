SUBROUTINE MACRO

USE SIMPLE
IMPLICIT NONE

    REAL SUM
    DO J=SP,NY
    DO I=0,NX
        SUM=0.0
        DO K=0,8
        SUM=SUM+F(K,I,J)
        END DO
        C(I,J)=SUM+DT*QXY(I,J)/2.0
    END DO
    END DO
    
    !C(:,0:SP-1)=0.0
RETURN
END