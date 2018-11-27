SUBROUTINE STREAM
    USE SIMPLE
    IMPLICIT NONE
    REAL FF(0:8,0:NX,0:NY)
    
    FF(:,:,:)=F(:,:,:)
    !!!!!!!!!!!!!! STREAMING!!!!!!!!!!!!!!!!!
    DO J=NY,0,-1
    DO I=NX,1,-1
    F(1,I,J)=FF(1,I-1,J)
    END DO
    END DO
    DO J=NY,1,-1
    DO I=1,NX
        F(2,I,J)=FF(2,I,J-1)
    END DO
    ENDDO
    DO J=0,NY
    DO I=0,NX-1
    F(3,I,J)=FF(3,I+1,J)
    END DO
    END DO
    DO J=0,NY-1
    DO I=NX,0,-1
        F(4,I,J)=FF(4,I,J+1)
    END DO
    END DO   
    DO J=NY,1,-1
        DO I=NX,1,-1
        F(5,I,J)=FF(5,I-1,J-1)
        END DO
    END DO
    DO J=NY,1,-1
        DO I=0,NX-1
        F(6,I,J)=FF(6,I+1,J-1)
        END DO
    END DO
    DO J=0,NY-1
        DO I=0,NX-1
        F(7,I,J)=FF(7,I+1,J+1)
        END DO
    END DO
    DO J=0,NY-1
        DO I=NX,2,-1
        F(8,I,J)=FF(8,I-1,J+1)
        END DO
    END DO

RETURN
END 