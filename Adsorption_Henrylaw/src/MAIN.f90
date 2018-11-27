!========================================================
!LBM -ADVECTION-DIFFUSION D2Q9

PROGRAM MAIN
USE SIMPLE
IMPLICIT NONE
REAL TIME

CALL INITIAL

STEP=0
TIME=0.0
OPEN(5,FILE='LOG.PLT',ACCESS='APPEND')
WRITE(5,*) "VARIABLES=TIME,ERROR,ERROR1,Qs_ave"
!================   MAIN LOOP  ====================!
DO WHILE(STEP*DT*UMAX/LY.LE.100)
    STEP=STEP+1
    TIME=STEP*DT*UMAX/LY
    PRINT *,TIME,' ',ERROR,ERROR1,QS(NS/2)
    WRITE(5,*)TIME,ERROR,ERROR1,SUM(QS)/NS
        
        DO I=0,NX
        DO J=SP,NY
            IF(J.EQ.SP)THEN 
                DCDN(I,J)=(3.0*C(I,SP)-4.0*C(I,SP+1)+C(I,SP+2))/(-2.0*DY)!(C(I,SP)-C(I,SP+1))/(-DY)!
                !DCDN([0,NX],J)=(C([0,NX],SP)-C([0,NX],SP+1))/(-DY)
                
            ELSEIF(J.LE.NY-1)THEN
                DCDN(I,J)=(C(I,J+1)-C(I,J-1))/(2.0*DY)!(3.0*C(I,J)-4.0*C(I,J+1)+C(I,J+2))/(2.0*DY)!
            ELSE
                DCDN(I,J)=(3.0*C(I,J)-4.0*C(I,J-1)+C(I,J-2))/(2.0*DY)
                !DCDN(I,J)=(C(I,J)-C(I,J-1))/DY
            ENDIF 
        ENDDO
        ENDDO
    
        CALL IB_ACTION
        
        CALL FEQ_COLLID
    
        CALL STREAM
        
        CALL BOUND_IB
        !CALL BOUND
        
        CALL MACRO           
        
        IF(MOD(STEP,1000).EQ.0)THEN
            OPEN(3,FILE='BOTTFLUX_LB.PLT')!
            WRITE(3,*)'VARIABLES=X/L,LB_RES,LEVESQUE'
            DO I=1,NX
            WRITE(3,*)X(I)/LX,CS_LB(I),LEVESQUE_LB(I)
            ENDDO
            CLOSE(3)
            OPEN(3,FILE='BOTTfLUX_ib.PLT')!
            WRITE(3,*)'VARIABLES=X/L,IB_RES,LEVESQUE'
            DO I=2,NS
            WRITE(3,*)XS(I)/LX,CS_IB(I),LEVESQUE(I)
            ENDDO
            CLOSE(3)
            
            OPEN(3,FILE='DCDN.PLT')!
            WRITE(3,*)'zone, T=DCDN'
            DO I=1,NX
            WRITE(3,*)X(I)/LX,DCDN(I,SP)
            ENDDO
            WRITE(3,*)'zone, T=DCSDN'
            DO I=1,NS
            WRITE(3,*)Xs(i)/lx, DCSDN(I)
            ENDDO
            CLOSE(3)
            
            CALL OUTPUT
            !PAUSE
        ENDIF
    
    
ENDDO
    
    
END
!======================================================