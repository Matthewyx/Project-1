! BOUNDARY CONDITIONS
    ! LEFT BOUNDARY CONDITION, THE TEMPERATURE IS GIVEN, C0
SUBROUTINE BOUND_IB
USE SIMPLE
IMPLICIT NONE
!EXTERNAL FUN_FEQ
!!!XIAOYI HE 2002
!REAL LABTA
!LABTA=KC/CS**2*(1.0+1.0/TAU)   (xiaoyi He. reactive transport)

    CALL ERROR_CALCULATOR1
    
    !LEFT BOUNDARY CONDITION
    DO J=1,NY
        F(1,0,J)=W(1)*C0+W(3)*C0-F(3,0,J)
        F(5,0,J)=W(5)*C0+W(7)*C0-F(7,0,J)
        F(8,0,J)=W(8)*C0+W(6)*C0-F(6,0,J)
    END DO
    !RIGHT BOUNDARY CONDITION
    DO J=0,NY
        F(3,NX,J)=F(3,NX-1,J)
        F(6,NX,J)=F(6,NX-1,J)
        F(7,NX,J)=F(7,NX-1,J)
    ENDDO
    !TOP BOUNDARY CONDITIONS
    DO I=0,NX
        !F(1,I,NY)=-F(3,I,NY)
        F(4,I,NY)=F(4,I,NY-1)!-F(2,I,NY)
        F(7,I,NY)=F(7,I,NY-1)!-F(5,I,NY)
        F(8,I,NY)=F(8,I,NY-1)!-F(6,I,NY)
        !F(0,I,NY)=0
        
        F(2,I,0)=-F(4,I,0)
        F(5,I,0)=-F(7,I,0)
        F(6,I,0)=-F(8,I,0)
        F(1,I,0)=-F(3,I,0)
        F(0,I,0)=0.0
    ENDDO
    !BOTTOM BOUNDARY CONDITIONS,KANG ETAL.
        

RETURN
END
    
SUBROUTINE IB_ACTION
USE SIMPLE
IMPLICIT NONE
REAL DELTX,DELTY

!!!!!!!!!!!!!!!!!!!!!!!!!!Interpolation for Lagarangian parameters from Eularian nodes!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO I=1,NS
    IF(XS(I).LE.X(1))THEN
      X1=0
      X2=3
    ELSEIF(XS(I).LE.X(3))THEN
        X1=(XS(I)-X(0))/DX+1-1
        X2=(XS(I)-X(0))/DX+1+2
    ELSEIF(XS(I).LE.X(NX-3))THEN
        X1=(XS(I)-X(0))/DX-3
        X2=(XS(I)-X(0))/DX+3
    ELSE
        X1=NX-3!(XS(I)-X(0))/DX-3!
        X2=NX
    ENDIF
        Y1=(YS(I)-Y(0))/DY+1-1
        Y2=(YS(I)-Y(0))/DY+1+2

        DCSDN(i)=0.0
        CIB(i)  =0.0
        DO M=X1,X2
        DO N=Y1,Y2
            RX=ABS(XS(I)-X(M))/DX
            RY=ABS(YS(I)-Y(N))/DY
            
            !!!!!3-POINT DIRAC DELTA FUNCTION!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !    IF(ABS(RX).LE.0.5)THEN
            !        DELTX=(1+SQRT(1-3.0*RX*RX))/3.0
            !    ELSEIF(ABS(RX).LE.1.5)THEN
            !        DELTX=(5-3.0*ABS(RX)-SQRT(-2.0+6*ABS(RX)-3.0*RX*RX))/6.0
            !    ELSE
            !        DELTX=0.0
            !    ENDIF
            !    
            !    IF(ABS(RY).LE.0.5)THEN
            !        DELTY=(1+SQRT(1-3.0*RY*RY))/3.0
            !    ELSEIF(ABS(RY).LE.1.5)THEN
            !        DELTY=(5-3.0*ABS(RY)-SQRT(-2.0+6*ABS(RY)-3.0*RY*RY))/6.0
            !    ELSE
            !        DELTY=0.0
            !   ENDIF
            !!!!!4-POINT DIRAC DELTA FUNCTION!!!!!!!!!!!!!!
            IF((RX.LE.2.0).AND.(RY.LE.2.0))THEN
                DELTX=(1.0+COS(0.5*PI*RX))/4.0
                DELTY=(1.0+COS(0.5*PI*RY))/4.0
                CIB(i)=CIB(i)+C(M,N)*DELTX*DELTY*1.35
                DCSDN(I)=DCSDN(I)+DCDN(M,N)*DELTX*DELTY*1.35
            ENDIF  
            
        ENDDO
        ENDDO
        
ENDDO  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Feedback law!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
QS(:)=0.0
DO I=1,NS
    QS(I)=2.0*(DS*DCSDN(I)-KC*CIB(I))
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Diffusing the feedback force to ambient Eularian points!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
QXY(:,:)=0.0
DO I=1,NS
    IF(XS(I).LE.X(1))THEN
      X1=0
      X2=3
    ELSEIF(XS(I).LE.X(3))THEN
        X1=(XS(I)-X(0))/DX+1-1
        X2=(XS(I)-X(0))/DX+1+2
    ELSEIF(XS(I).LE.X(NX-3))THEN
        X1=(XS(I)-X(0))/DX-3
        X2=(XS(I)-X(0))/DX+3
    ELSE
        X1=NX-3!(XS(I)-X(0))/DX-3!
        X2=NX
    ENDIF
        Y1=(YS(I)-Y(0))/DY+1-1
        Y2=(YS(I)-Y(0))/DY+1+2

        DO M=X1,X2
        DO N=Y1,Y2
            RX= ABS(XS(I)-X(M))/DX
            RY= ABS(YS(I)-Y(N))/DY
            !!!!!!!3-POINT DIRAC DELTA FUNCTION!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !    IF(ABS(RX).LE.0.5)THEN
            !        DELTX=(1+SQRT(1-3.0*RX*RX))/3.0
            !    ELSEIF(ABS(RX).LE.1.5)THEN
            !        DELTX=(5-3.0*ABS(RX)-SQRT(-2.0+6*ABS(RX)-3.0*RX*RX))/6.0
            !    ELSE
            !       DELTX=0.0
            !    ENDIF
            !
            !    IF(ABS(RY).LE.0.5)THEN
            !        DELTY=(1+SQRT(1-3.0*RY*RY))/3.0
            !    ELSEIF(ABS(RY).LE.1.5)THEN
            !        DELTY=(5-3.0*ABS(RY)-SQRT(-2.0+6*ABS(RY)-3.0*RY*RY))/6.0
            !    ELSE
            !        DELTY=0.0
            !    ENDIF           
            !!!!!4-POINT DIRAC DELTA FUNCTION!!!!!!!!!!!!!!
                IF((RX.LE.2.0).AND.(RY.LE.2.0))THEN
                    DELTX=(1.0+COS(0.5*PI*RX))/4.0
                    DELTY=(1.0+COS(0.5*PI*RY))/4.0
                    QXY(M,N)=QS(I)*ARC_L/DX**2*DELTX*DELTY  
                ENDIF
                 
        ENDDO
        ENDDO
ENDDO
    
RETURN
END
    
SUBROUTINE ERROR_CALCULATOR1
USE SIMPLE

IMPLICIT NONE

REAL TEMP1,TEMP2
!!!!!!!!!!!!!ERROR CHECK PROCESS!!!!!!!!!!!!!
        DO I=2,NS
            CS_IB(I)=((LY-SP*DY)/C0)*DCSDN(I)
            
            TEMP1=UMAX*(LY-SP*DY)**2
            TEMP2=XS(I)*DS
            LEVESQUE(I)=0.854*(TEMP1/TEMP2)**(1.0/3.0)
        
            ERROR=ERROR+ABS(CS_IB(I)-LEVESQUE(I))!**2
        END DO
        ERROR=ERROR/SUM(LEVESQUE(2:NS))
        !ERROR=SQRT(ERROR/(NS-1))/(SUM(LEVESQUE(2:NS))/(NS-1))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO I=1,NX
            CS_LB(I)=((LY-SP*DY)/C0)*DCDN(I,SP)
            TEMP1=(2.0*UMAX/3.0)*(LY-SP*DY)**2
            TEMP2=X(I)*DS
            LEVESQUE_LB(I)=0.854*(TEMP1/TEMP2)**(1.0/3.0)
            
            ERROR1=ERROR1+ABS(CS_LB(I)-LEVESQUE_LB(I))!**2
        ENDDO
        ERROR1=ERROR1/SUM(LEVESQUE_LB(1:NX))
RETURN
END 

!!!!INTERPOLATION OF THE REQUIRED VALUE ON LAGARANGIAN NODES!!!!!!!!!!!!!!
!!!!FROM FLUID TO VIRTUAL SOLID NODES/LAGARANGS
SUBROUTINE DELT_EULA2LAG(XX)
USE SIMPLE
INTEGER XX
REAL DELTX,DELTY

DO I=1,NS
    IF(XS(I).LE.X(1))THEN
      X1=0
      X2=3
    ELSEIF(XS(I).LE.X(NX-2))THEN
        X1=(XS(I)-X(0))/DX+1-1
        X2=(XS(I)-X(0))/DX+1+2
    ELSE
        X1=NX-3
        X2=NX
    ENDIF
        Y1=(YS(I)-Y(0))/DY-1
        Y2=(YS(I)-Y(0))/DY+2

        DCSDN(i)=0.0
        CIB(i)  =0.0
        DO M=X1,X2
        DO N=Y1,Y2
            RX=(XS(I)-X(M))/DX!ABS(XS(I)-X(M))/DX
            RY=(YS(I)-Y(N))/DY!ABS(YS(I)-Y(N))/DY
            IF(XX.EQ.3)THEN
            !!!!!3-POINT DIRAC DELTA FUNCTION!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                IF(ABS(RX).LE.0.5)THEN
                    DELTX=(1+SQRT(1-3.0*RX*RX))/3.0
                ELSEIF(ABS(RX).LE.1.5)THEN
                    DELTX=(5-3.0*ABS(RX)-SQRT(-2.0+6*ABS(RX)-3.0*RX*RX))/6.0
                ELSE
                    DELTX=0.0
                ENDIF
                
                IF(ABS(RY).LE.0.5)THEN
                    DELTY=(1+SQRT(1-3.0*RY*RY))/3.0
                ELSEIF(ABS(RY).LE.1.5)THEN
                    DELTY=(5-3.0*ABS(RY)-SQRT(-2.0+6*ABS(RY)-3.0*RY*RY))/6.0
                ELSE
                    DELTY=0.0
                ENDIF
            ELSEIF(XX.EQ.4)THEN
            !!!!!4-POINT DIRAC DELTA FUNCTION!!!!!!!!!!!!!!
                IF((RX.LE.2.0).AND.(RY.LE.2.0))THEN
                    DELTX=(1.0+COS(0.5*PI*RX))/4.0
                    DELTY=(1.0+COS(0.5*PI*RY))/4.0
                ENDIF  
            ENDIF
            
        CIB(i)=CIB(i)+C(M,N)*DELTX*DELTY
        DCSDN(I)=DCSDN(I)+DCDN(M,N)*DELTX*DELTY
        ENDDO
        ENDDO
ENDDO  

END
    
    
SUBROUTINE DELT_LAG2EULA(YY)
USE SIMPLE
INTEGER YY
REAL DELTX,DELTY

QXY(:,:)=0.0
DO I=1,NS
    IF(XS(I).LE.X(1))THEN
      X1=0
      X2=3
    ELSEIF(XS(I).LE.X(NX-2))THEN
        X1=(XS(I)-X(0))/DX+1-1
        X2=(XS(I)-X(0))/DX+1+2
    ELSE
        X1=NX-3
        X2=NX
    ENDIF
        Y1=(YS(I)-Y(0))/DY-1
        Y2=(YS(I)-Y(0))/DY+2
        
        DO M=X1,X2
        DO N=Y1,Y2
            RX=(XS(I)-X(M))/DX!ABS(XS(I)-X(M))/DX
            RY=(YS(I)-Y(N))/DY!ABS(YS(I)-Y(N))/DY
            
            IF(YY.EQ.3)THEN
                !!!!!3-POINT DIRAC DELTA FUNCTION!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                IF(ABS(RX).LE.0.5)THEN
                    DELTX=(1+SQRT(1-3.0*RX*RX))/3.0
                ELSEIF(ABS(RX).LE.1.5)THEN
                    DELTX=(5-3.0*ABS(RX)-SQRT(-2.0+6*ABS(RX)-3.0*RX*RX))/6.0
                ELSE
                    DELTX=0.0
                ENDIF
            
                IF(ABS(RY).LE.0.5)THEN
                    DELTY=(1+SQRT(1-3.0*RY*RY))/3.0
                ELSEIF(ABS(RY).LE.1.5)THEN
                    DELTY=(5-3.0*ABS(RY)-SQRT(-2.0+6*ABS(RY)-3.0*RY*RY))/6.0
                ELSE
                    DELTY=0.0
                ENDIF           
            ELSEIF(YY.EQ.4)THEN
                !!!!!4-POINT DIRAC DELTA FUNCTION!!!!!!!!!!!!!!
                IF((RX.LE.2.0).AND.(RY.LE.2.0))THEN
                    DELTX=(1.0+COS(0.5*PI*RX))/4.0
                    DELTY=(1.0+COS(0.5*PI*RY))/4.0
                ELSE
                    DELTX=0.0
                    DELTY=0.0
                ENDIF
            ENDIF
            
            QXY(M,N)=QS(I)*ARC_L/DX**2*DELTX*DELTY       
        ENDDO
        ENDDO
ENDDO

    END
