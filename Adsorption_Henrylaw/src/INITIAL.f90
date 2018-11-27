SUBROUTINE INITIAL

USE SIMPLE

IMPLICIT NONE
REAL,EXTERNAL:: FUN_FEQ

PI=4.0*ATAN(1.0)
!!!!!!!!!!!!!!!!!!!!
DX=LX/NX
DY=Ly/Ny
DT=DX
!!!!!!!!!!!!!!!!!!MESH GRID!!!!!!!!!!!!!!!!
DO I=0,NX
    X(I)=I*DX
END DO
DO J=0,NY
    Y(J)=J*DY
END DO
!!!!!!!!!!VELOCITY DISTRIBUTION!!!!!!!!!!!
UMAX=0.06
DO I=0,NX
DO J=SP,NY
    U(I,J)=-4.0*UMAX*(Y(J)-SP*DY)*(Y(J)-LY)/(LY-SP*dy)**2
ENDDO
ENDDO
V=0.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CS=SQRT(1.0/3.0)*(DX/DT)
KC=1.0          !HENRY'S CONSTANT
C0=1.0
DS=1.0/6.0
CK=DX/DT
CSQ=CK*CK
TAU=3.0*DS/DT+0.5
OMEGA=1./TAU

W(0)=4./9.
CX(0)=0.0
CY(0)=0.0
DO I=1,4
CX(I)=COS((I-1)*0.5*PI)*DX/DT
CY(I)=SIN((I-1)*0.5*PI)*DY/DT
W(I)=1./9.
ENDDO
DO I=5,8
CX(I)=SQRT(2.0)*COS(0.25*PI+0.5*PI*(I-5))*DX/DT
CY(I)=SQRT(2.0)*SIN(0.25*PI+0.5*PI*(I-5))*DY/DT
W(I)=1./36.
ENDDO

C(:,:)=0.0
C(0,SP:NY)=C0
DO J=SP,NY
DO I=0,NX
    DO K=0,8
    F(K,I,J)=FUN_FEQ(K,U(I,J),V(I,J),C(I,J))
    ENDDO
ENDDO
ENDDO

!!!!!!!!!!!!!LAGARANGIAN NODES!!!!!!!!!!!!!!!!!!!

ARC_L=LX/(NS-1)
DO I=1,NS
    XS(I)=(I-1)*ARC_L
    YS(I)=SP*DY
ENDDO

QS =0.0
QXY =0.0
    
RETURN
END