PROGRAM MAIN
PARAMETER(IM=2401,JM=801,KM=9)
PARAMETER(ISY=400,ISB=320)
!INTEGER IM,JM,KM, I,J
INTEGER NSTEP
DOUBLE PRECISION ERROR, RE, DEN, U0
DOUBLE PRECISION X,Y
DOUBLE PRECISION DT, PI, TAO, CS
DOUBLE PRECISION UC, VC
DOUBLE PRECISION U,V,D
DOUBLE PRECISION FEQ,F
DOUBLE PRECISION CD,CL
INTEGER TIME(8)
CHARACTER (LEN=10) BIG_BEN(3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE PRECISION SU,SV,SX,SY,CENU,CENV,CENX,CENY,SFX,SFY,SPRKY,SXH,SYH
DOUBLE PRECISION FX,FY
DOUBLE PRECISION SUB,SVB,SXB,SYB,SFXB,SFYB,SXHB,SYHB
DOUBLE PRECISION DS
COMMON/DS/DS
COMMON/IBM1/SUB(ISB),SVB(ISB),SXB(ISB),SYB(ISB),SFXB(ISB),SFYB(ISB),SXHB(ISB),SYHB(ISB)
COMMON/FORCE/FX(IM,JM),FY(IM,JM)
COMMON/IBM/SU(ISY),SV(ISY),SX(ISY),SY(ISY),CENU,CENV,CENX,CENY,SFX(ISY),SFY(ISY),SPRKY,SXH(ISY),SYH(ISY)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
COMMON/STEP/NSTEP 
COMMON/PHT/FEQ(IM,JM,KM),F(IM,JM,KM)
COMMON/UVD/U(IM,JM),V(IM,JM),D(IM,JM)
COMMON/CONST1/DT, PI, TAO, CS
COMMON/CONST2/UC(KM), VC(KM)
COMMON/CONST/ERROR, RE, DEN, U0
COMMON/XY/X(IM,JM),Y(IM,JM)



OPEN(1, FILE='INPUT.DAT',STATUS='OLD')
READ(1,*)ERROR, RE, DEN, U0, SPRKY

DO I=1,IM
DO J=1,JM
 X(I,J)=30.0*(I-1)/((IM-1)*1.0)
 Y(I,J)=10.0*(J-1)/((JM-1)*1.0)
ENDDO
ENDDO

 
 WRITE(*,*)"PLEASE INPUT A NUMBER, 1 TO START, 2 TO CONTINUE!"
 READ(1,*) N
CLOSE(1)
 IF(N.EQ.1)THEN
 		NSTEP=0
        CALL INITIAL
 
   
   100  NSTEP=NSTEP+1
   
        IF(MOD(NSTEP,1000).EQ.0)THEN
        OPEN(1,FILE='TIME.TXT',ACCESS='APPEND')
 	    CALL DATE_AND_TIME(BIG_BEN(1),BIG_BEN(2),BIG_BEN(3),TIME)
        WRITE(1,*)NSTEP,'   ', BIG_BEN(1), BIG_BEN(2)
        CLOSE(1)
        ELSE
        ENDIF

		
      
		CALL ACTION1
		CALL COLLID
		CALL STREAM
		CALL MACRO
  		CALL PINGHT

        IF(MOD(NSTEP,10000).EQ.0)THEN
         CALL OUTPUT
		
        ELSE
        ENDIF

		IF(MOD(NSTEP,30000).EQ.0)THEN
          OPEN(1,FILE='F.TXT')
 	      DO I=1,IM
		  DO J=1,JM
		  DO K=1,KM
		    WRITE(1,*) F(I,J,K), FEQ(I,J,K)
		  ENDDO
		  ENDDO
		  ENDDO
          CLOSE(1)

		  OPEN(1,FILE='UVD.TXT')
 	      DO I=1,IM
		  DO J=1,JM
		  	WRITE(1,*) U(I,J), V(I,J), D(I,J)
		  ENDDO
		  ENDDO
          CLOSE(1)


		  OPEN(1,FILE='CONST.TXT')
 	      WRITE(1,*)DT, PI, TAO, CS
		  WRITE(1,*)NSTEP
		  DO I=1,KM
		  WRITE(1,*)UC(I),VC(I)
		  ENDDO
		  DO I=1,ISY
		  WRITE(1,*)SX(I),SY(I),SXH(I),SYH(I)
		  ENDDO
		  WRITE(1,*)CENX,CENY
		  DO I=1,ISB
		  WRITE(1,*)SXB(I),SYB(I),SXHB(I),SYHB(I)
		  ENDDO

		  CLOSE(1)
 
		
        ELSE
		ENDIF

        !!!!!!!!!!!!!CdCl
		!HERE WE USE THE CYLINDER DIAMETER
		
		  CD=0.0
          CL=0.0
		  DO I=1,ISY
           CD=CD-2.0*SFX(I)*DS/(1.0*U0**2.0*1.0)
          ENDDO
		  DO I=1,ISY
           CL=CL-2.0*SFY(I)*DS/(1.0*U0**2.0*1.0)
          ENDDO
		  
		  OPEN(1,FILE='CD.TXT',ACCESS='APPEND')
          WRITE(1,*)NSTEP,CD
		  WRITE(*,*)NSTEP,CD,CL
          CLOSE(1)

		  OPEN(1,FILE='CL.TXT',ACCESS='APPEND')
          WRITE(1,*)NSTEP,CL
		  CLOSE(1)


		  OPEN(1,FILE='CEN.TXT',ACCESS='APPEND')
          WRITE(1,*)NSTEP,CENX, CENY
          CLOSE(1)

	    
   		GOTO 100 
		
	ELSE

		  OPEN(1,FILE='F.TXT')
 	      DO I=1,IM
		  DO J=1,JM
		  DO K=1,KM
		    READ(1,*) F(I,J,K), FEQ(I,J,K)
		  ENDDO
		  ENDDO
		  ENDDO
          CLOSE(1)

		  OPEN(1,FILE='UVD.TXT')
 	      DO I=1,IM
		  DO J=1,JM
		  	READ(1,*) U(I,J), V(I,J), D(I,J)
		  ENDDO
		  ENDDO
          CLOSE(1)


		  OPEN(1,FILE='CONST.TXT')
 	      READ(1,*)DT, PI, TAO, CS
		  READ(1,*)NSTEP
		  DO I=1,KM
		  READ(1,*)UC(I),VC(I)
		  ENDDO
		  DO I=1,ISY
		  READ(1,*)SX(I),SY(I),SXH(I),SYH(I)
		  ENDDO
		  READ(1,*)CENX,CENY
		  DO I=1,ISB
		  READ(1,*)SXB(I),SYB(I),SXHB(I),SYHB(I)
		  ENDDO
		  
		  CLOSE(1)

 200  NSTEP=NSTEP+1
   
        IF(MOD(NSTEP,1000).EQ.0)THEN
        OPEN(1,FILE='TIME.TXT',ACCESS='APPEND')
 	    CALL DATE_AND_TIME(BIG_BEN(1),BIG_BEN(2),BIG_BEN(3),TIME)
        WRITE(1,*) NSTEP, '   ', BIG_BEN(1), BIG_BEN(2)
        CLOSE(1)
        ELSE
        ENDIF

		
      
		CALL ACTION1
		CALL COLLID
		CALL STREAM
		CALL MACRO
  		CALL PINGHT

        IF(MOD(NSTEP,200).EQ.0)THEN
         CALL OUTPUT
		ELSE
        ENDIF

		!!!!!!!!!!!!!CdCl
		!HERE WE USE THE CYLINDER DIAMETER
		
		  CD=0.0
          CL=0.0
		  DO I=1,ISY
           CD=CD-2.0*SFX(I)*dt**2.0/(1.0*U0**2.0*1.0)
          ENDDO
		  DO I=1,ISY
           CL=CL-2.0*SFY(I)*dt**2.0/(1.0*U0**2.0*1.0)
          ENDDO
		  
		  OPEN(1,FILE='CD.TXT',ACCESS='APPEND')
          WRITE(1,*)NSTEP,CD
		  WRITE(*,*)NSTEP,CD,CL
          CLOSE(1)

		  OPEN(1,FILE='CL.TXT',ACCESS='APPEND')
          WRITE(1,*)NSTEP,CL
		  CLOSE(1)


		  OPEN(1,FILE='CEN.TXT',ACCESS='APPEND')
          WRITE(1,*)NSTEP,CENX, CENY
          CLOSE(1)

	    
   		GOTO 200 
 


		ENDIF
	
 END

SUBROUTINE INITIAL
PARAMETER(IM=2401,JM=801,KM=9)
PARAMETER(ISY=400,ISB=320)
!INTEGER IM,JM,KM,I,J,K
DOUBLE PRECISION ERROR, RE, DEN, U0
DOUBLE PRECISION DT, PI, TAO, CS
DOUBLE PRECISION UC, VC
DOUBLE PRECISION U,V,D
DOUBLE PRECISION FEQ,F
DOUBLE PRECISION X,Y
!!!!!!!!!!!!!!!!!!
DOUBLE PRECISION SU,SV,SX,SY,CENU,CENV,CENX,CENY,SFX,SFY,SPRKY,SXH,SYH
DOUBLE PRECISION FX,FY
DOUBLE PRECISION SUB,SVB,SXB,SYB,SFXB,SFYB,SXHB,SYHB
DOUBLE PRECISION DS
COMMON/DS/DS
COMMON/IBM1/SUB(ISB),SVB(ISB),SXB(ISB),SYB(ISB),SFXB(ISB),SFYB(ISB),SXHB(ISB),SYHB(ISB)
COMMON/FORCE/FX(IM,JM),FY(IM,JM)
COMMON/IBM/SU(ISY),SV(ISY),SX(ISY),SY(ISY),CENU,CENV,CENX,CENY,SFX(ISY),SFY(ISY),SPRKY,SXH(ISY),SYH(ISY)
!!!!!!!!!!!!!!!!!!
COMMON/XY/X(IM,JM),Y(IM,JM)
COMMON/PHT/FEQ(IM,JM,KM),F(IM,JM,KM)
COMMON/UVD/U(IM,JM),V(IM,JM),D(IM,JM)
COMMON/CONST1/DT, PI, TAO, CS
COMMON/CONST2/UC(KM), VC(KM)
COMMON/CONST/ERROR, RE, DEN, U0

PI=4.0*ATAN(1.0)
DT=30.0/(IM-1)
CS=1.0/SQRT(3.0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TAO=0.5*(1.0+6.0*U0*1.0/(RE*DT))
WRITE(*,*)'TAO IS:',TAO


DO I=1,KM
  IF(I.LE.4)THEN
  UC(I)=COS((I-1)*0.5*PI)
  VC(I)=SIN((I-1)*0.5*PI)
  ELSE IF(I.LE.8)THEN
  UC(I)=SQRT(2.0)*COS(0.25*PI+0.5*PI*(I-5))
  VC(I)=SQRT(2.0)*SIN(0.25*PI+0.5*PI*(I-5))
  ELSE
  UC(I)=0.0
  VC(I)=0.0
  ENDIF
ENDDO

DO I=1,IM
DO J=1,JM
 U(I,J)=U0
 V(I,J)=0
 D(I,J)=DEN
ENDDO
ENDDO

CALL PINGHT

DO I=1,IM
DO J=1,JM
DO K=1,KM
 F(I,J,K)=FEQ(I,J,K)
ENDDO
ENDDO
ENDDO

!!!!!!!!!!LAGRANGE!!!!!!!!!!!

DO I=1,ISB
SXB(I)=8.0+2.5*(I-1)/(ISB-1)
SXHB(I)=SXB(I)
SYB(I)=5.0125
SYHB(I)=SYB(I)
SUB(I)=0.0
SVB(I)=0.0
SFXB(I)=0.0
SFYB(I)=0.0
ENDDO

DO I=1,ISY
SX(I)=5.0+0.5*COS(2*PI*(I-1)*1.0/(ISY*1.0))
SXH(I)=SX(I)
SY(I)=5.0+0.5*SIN(2*PI*(I-1)*1.0/(ISY*1.0))
SYH(I)=SY(I)
SU(I)=0.0
SV(I)=0.0
SFX(I)=0.0
SFY(I)=0.0
ENDDO
DS=Pi*2.0*0.5/ISY

DO I=1,IM
DO J=1,JM
FX(I,J)=0.0
FY(I,J)=0.0
ENDDO
ENDDO
CENU=0
CENV=0
CENX=5.0
CENY=5.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
RETURN
END

SUBROUTINE PINGHT
PARAMETER(IM=2401,JM=801,KM=9)
INTEGER IM,JM,KM,I,J,K
DOUBLE PRECISION ERROR, RE, DEN, U0
DOUBLE PRECISION DT, PI, TAO, CS
DOUBLE PRECISION UC, VC
DOUBLE PRECISION U,V,D
DOUBLE PRECISION FEQ,F
COMMON/PHT/FEQ(IM,JM,KM),F(IM,JM,KM)
COMMON/UVD/U(IM,JM),V(IM,JM),D(IM,JM)
COMMON/CONST1/DT, PI, TAO, CS
COMMON/CONST2/UC(KM), VC(KM)
COMMON/CONST/ERROR, RE, DEN, U0

  DO I=1,IM
  DO J=1,JM
    DO K=1, KM
	IF(K.LE.4)THEN
	FEQ(I,J,K)=(D(I,J)+3.0*(U(I,J)*UC(K)+V(I,J)*VC(K))+4.5*(U(I,J)*UC(K)+V(I,J)*VC(K))**2.0-1.5*(U(I,J)**2.0+V(I,J)**2.0))/9.0
	ELSE IF(K.LE.8)THEN
	FEQ(I,J,K)=(D(I,J)+3.0*(U(I,J)*UC(K)+V(I,J)*VC(K))+4.5*(U(I,J)*UC(K)+V(I,J)*VC(K))**2.0-1.5*(U(I,J)**2.0+V(I,J)**2.0))/36.0
	ELSE
	FEQ(I,J,K)=(D(I,J)-1.5*(U(I,J)**2.0+V(I,J)**2.0))*4.0/9.0
	ENDIF
	ENDDO
  ENDDO
  ENDDO
  RETURN
END


SUBROUTINE ACTION1
PARAMETER(IM=2401,JM=801,KM=9)
PARAMETER(ISY=400,ISB=320)
INTEGER IM,JM,KM, I, J, I1,I2,J1,J2, M, N
DOUBLE PRECISION ERROR, RE, DEN, U0
DOUBLE PRECISION X,Y
DOUBLE PRECISION DT, PI, TAO, CS
DOUBLE PRECISION U,V,D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE PRECISION SU,SV,SX,SY,CENU,CENV,CENX,CENY,SFX,SFY,SPRKY,SXH,SYH
DOUBLE PRECISION FX,FY,SPRKY_a,SPRKY_b
!DOUBLE PRECISION SUB,SVB,SXB,SYB,SFXB,SFYB,SXHB,SYHB,SPRKB,SPRKBB
!DOUBLE PRECISION LEN0, LENDEF(ISB)
!COMMON/IBM1/SUB(ISB),SVB(ISB),SXB(ISB),SYB(ISB),SFXB(ISB),SFYB(ISB),SXHB(ISB),SYHB(ISB)
DOUBLE PRECISION DS
COMMON/DS/DS
COMMON/FORCE/FX(IM,JM),FY(IM,JM)
COMMON/IBM/SU(ISY),SV(ISY),SX(ISY),SY(ISY),CENU,CENV,CENX,CENY,SFX(ISY),SFY(ISY),SPRKY,SXH(ISY),SYH(ISY)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
COMMON/UVD/U(IM,JM),V(IM,JM),D(IM,JM)
COMMON/CONST1/DT, PI, TAO, CS
COMMON/CONST/ERROR, RE, DEN, U0
COMMON/XY/X(IM,JM),Y(IM,JM)
COMMON/STEP/NSTEP

SPRKY_a=0.0*SPRKY*U0**2 !50.0*SPRKY*U0**2

SPRKY_b=SPRKY*U0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CYNLINDER

DO I=1,ISY
I1=SXH(I)/DT-3
I2=SXH(I)/DT+3
J1=SYH(I)/DT-3
J2=SYH(I)/DT+3
  SU(I)=0.0
  SV(I)=0.0
  DO M=I1,I2
  DO N=J1,J2
  IF(ABS(SXH(I)-X(M,N))/DT.LE.2.0.AND.ABS(SYH(I)-Y(M,N))/DT.LE.2.0)THEN
  SU(I)=SU(I)+U(M,N)*(1.0+COS(0.5*PI*ABS(SXH(I)-X(M,N))/DT))*(1.0+COS(0.5*PI*ABS(SYH(I)-Y(M,N))/DT))/16.0 
  SV(I)=SV(I)+V(M,N)*(1.0+COS(0.5*PI*ABS(SXH(I)-X(M,N))/DT))*(1.0+COS(0.5*PI*ABS(SYH(I)-Y(M,N))/DT))/16.0 
  ELSE
  ENDIF
  ENDDO
  ENDDO 
ENDDO

CENU=0
CENV=0
DO I=1,ISY
SXH(I)=SXH(I) !+DT*SU(I)
SYH(I)=SYH(I) !+DT*SV(I)

SFX(I)=SPRKY_a*(SX(I)-SXH(I))+SPRKY_b*(0.0-SU(I))
SFY(I)=SPRKY_a*(SY(I)-SYH(I))+SPRKY_b*(0.0-SV(I))
CENU=CENU+SU(I)
CENV=CENV+SV(I)
ENDDO
CENU=CENU/ISY
CENV=CENV/ISY
CENX=CENX+DT*CENU
CENY=CENY+DT*CENV

DO I=1,IM
DO J=1,JM
 FX(I,J)=0.0
 FY(I,J)=0.0
ENDDO
ENDDO 

DO I=1,ISY

I1=SXH(I)/DT-3
I2=SXH(I)/DT+3
J1=SYH(I)/DT-3
J2=SYH(I)/DT+3

  DO M=I1,I2
  DO N=J1,J2
  IF(ABS(SXH(I)-X(M,N))/DT.LE.2.0.AND.ABS(SYH(I)-Y(M,N))/DT.LE.2.0)THEN
  FX(M,N)=FX(M,N)+SFX(I)*(1.0+COS(0.5*PI*ABS(SXH(I)-X(M,N))/DT))*(1.0+COS(0.5*PI*ABS(SYH(I)-Y(M,N))/DT))/16.0 *DS/DT**2
  FY(M,N)=FY(M,N)+SFY(I)*(1.0+COS(0.5*PI*ABS(SXH(I)-X(M,N))/DT))*(1.0+COS(0.5*PI*ABS(SYH(I)-Y(M,N))/DT))/16.0 *DS/DT**2
  ELSE
  ENDIF
  ENDDO
  ENDDO 

ENDDO


END

SUBROUTINE STREAM
PARAMETER(IM=2401,JM=801,KM=9)
INTEGER IM,JM,KM, I,J,K
DOUBLE PRECISION ERROR, RE, DEN, U0
DOUBLE PRECISION DT, PI, TAO, CS
DOUBLE PRECISION UC, VC
DOUBLE PRECISION U,V,D
DOUBLE PRECISION FEQ,F
DOUBLE PRECISION FF(IM,JM,KM)
COMMON/PHT/FEQ(IM,JM,KM),F(IM,JM,KM)
COMMON/UVD/U(IM,JM),V(IM,JM),D(IM,JM)
COMMON/CONST1/DT, PI, TAO, CS
COMMON/CONST2/UC(KM), VC(KM)
COMMON/CONST/ERROR, RE, DEN, U0

   DO I=1,IM
   DO J=1,JM
   DO K=1,KM
    FF(I,J,K)=F(I,J,K)
   ENDDO
   ENDDO
   ENDDO

   DO J=1,JM
   DO I=2,IM
    F(I,J,1)=FF(I-1,J,1)
   ENDDO
   ENDDO

   DO J=2,JM
   DO I=1,IM
        F(I,J,2)=FF(I,J-1,2)
   ENDDO
   ENDDO

   DO J=1,JM
   DO I=1,IM-1
   
    F(I,J,3)=FF(I+1,J,3)
	
   ENDDO
   ENDDO

   DO J=1,JM-1
   DO I=1,IM
    
    F(I,J,4)=FF(I,J+1,4)
	
   ENDDO
   ENDDO

   DO J=2,JM
   DO I=2,IM
   
    F(I,J,5)=FF(I-1,J-1,5)
	
   ENDDO
   ENDDO

   DO J=2,JM
   DO I=1,IM-1
    
    F(I,J,6)=FF(I+1,J-1,6)
	
   ENDDO
   ENDDO

   DO J=1,JM-1
   DO I=1,IM-1
    
    F(I,J,7)=FF(I+1,J+1,7)
	
   ENDDO
   ENDDO

   DO J=1,JM-1
   DO I=2,IM
    
    F(I,J,8)=FF(I-1,J+1,8)
	
   ENDDO
   ENDDO

!BOUNDARY CONDITION
 
   DO K=1,KM
   DO I=1, IM
   F(I,1,K)=FEQ(I,1,K)
   F(I,JM,K)=FEQ(I,JM,K)
   ENDDO
   DO J=1, JM
   F(1,J,K)=FEQ(1,J,K)
   F(IM,J,K)=FEQ(IM,J,K)
   ENDDO
  ENDDO  
   
  ! INLET	
  !  DO J=2, JM-1
  !  F(1,J,1)=F(1,J,3)+2.0*U0/3.0
!	F(1,J,5)=F(1,J,7)-0.5*(F(1,J,2)-F(1,J,4))+U0/6.0
!	F(1,J,8)=F(1,J,6)+0.5*(F(1,J,2)-F(1,J,4))+U0/6.0
!	ENDDO

 ! upper and lower (periodic) 
   !DO K=1,KM
   !DO I=1,IM
   !F(I,1,K)=F(I,JM-1,K)
   !F(I,JM,K)=F(I,2,K)
   !ENDDO
 ! outlet
  ! DO J=1, JM
  ! F(IM,J,K)=2*F(IM-1,J,K)-F(IM-2,J,K)
  ! ENDDO
  ! ENDDO 
   
   ! INLET SINGULAR 
!	F(1,1,1)=F(1,1,3)
!	F(1,1,2)=F(1,1,4)
!	F(1,1,5)=F(1,1,7)
!	F(1,1,6)=0.5*(D(1,2)-F(1,1,9)-F(1,1,1)-F(1,1,2)-F(1,1,3)-F(1,1,4)-F(1,1,5)-F(1,1,7))
!	F(1,1,8)=0.5*(D(1,2)-F(1,1,9)-F(1,1,1)-F(1,1,2)-F(1,1,3)-F(1,1,4)-F(1,1,5)-F(1,1,7))
	

!	F(1,JM,1)=F(1,JM,3)
!	F(1,JM,4)=F(1,JM,2)
!	F(1,JM,8)=F(1,JM,6)
!	F(1,JM,5)=0.5*(D(1,JM-1)-F(1,JM,9)-F(1,JM,1)-F(1,JM,2)-F(1,JM,3)-F(1,JM,4)-F(1,JM,6)-F(1,JM,8))
!	F(1,JM,7)=0.5*(D(1,JM-1)-F(1,JM,9)-F(1,JM,1)-F(1,JM,2)-F(1,JM,3)-F(1,JM,4)-F(1,JM,6)-F(1,JM,8))
	  
	RETURN
END


SUBROUTINE MACRO
PARAMETER(IM=2401,JM=801,KM=9)
INTEGER IM,JM,KM, I,J,K
DOUBLE PRECISION ERROR, RE, DEN, U0
DOUBLE PRECISION UC, VC
DOUBLE PRECISION U,V,D
DOUBLE PRECISION UU,VV,DD
DOUBLE PRECISION FEQ,F
DOUBLE PRECISION FX,FY
DOUBLE PRECISION DT, PI, TAO, CS
COMMON/CONST1/DT, PI, TAO, CS
COMMON/FORCE/FX(IM,JM),FY(IM,JM)
COMMON/PHT/FEQ(IM,JM,KM),F(IM,JM,KM)
COMMON/UVD/U(IM,JM),V(IM,JM),D(IM,JM)
COMMON/CONST2/UC(KM), VC(KM)
COMMON/CONST/ERROR, RE, DEN, U0

  DO I=1,IM
  DO J=1,JM
   UU=0.0
   VV=0.0
   DD=0.0
   DO K=1,KM
   DD=DD+F(I,J,K)
   UU=UU+F(I,J,K)*UC(K)
   VV=VV+F(I,J,K)*VC(K)
   ENDDO

   D(I,J)=DD
   U(I,J)=(UU+0.5*DT*FX(I,J))/DEN
   V(I,J)=(VV+0.5*DT*FY(I,J))/DEN
  
  ENDDO
  ENDDO

  DO I=1,IM
   U(I,JM)=U0
   V(I,JM)=0.0
   D(I,JM)=DEN
   U(I,1)=U0
   V(I,1)=0.0
   D(I,1)=DEN
  ENDDO

  DO J=2,JM-1
   U(1,J)=U0
   V(1,J)=0.0
   D(1,J)=DEN

   U(IM,J)=U0
   V(IM,J)=0.0
   D(IM,J)=DEN
  ENDDO
  RETURN
END


SUBROUTINE COLLID
PARAMETER(IM=2401,JM=801,KM=9)
INTEGER IM,JM,KM, I,J,K
DOUBLE PRECISION ERROR, RE, DEN, U0
DOUBLE PRECISION DT, PI, TAO, CS
DOUBLE PRECISION FEQ,F
DOUBLE PRECISION FX,FY,U,V,D,UC,VC
COMMON/CONST2/UC(KM), VC(KM)
COMMON/UVD/U(IM,JM),V(IM,JM),D(IM,JM)
COMMON/FORCE/FX(IM,JM),FY(IM,JM)
COMMON/PHT/FEQ(IM,JM,KM),F(IM,JM,KM)
COMMON/CONST1/DT, PI, TAO, CS
COMMON/CONST/ERROR, RE, DEN, U0
 
 DO I=1,IM
 DO J=1,JM
 DO K=1,KM
   IF(K.LE.4)THEN
   F(I,J,K)=F(I,J,K)*(1.0-1.0/TAO)+FEQ(I,J,K)/TAO+DT*(1.0-0.5/TAO)*((3.0*(UC(K)-U(I,J))+9.0*(UC(K)*U(I,J)+VC(K)*V(I,J))*UC(K))*FX(I,J)+(3.0*(VC(K)-V(I,J))+9.0*(UC(K)*U(I,J)+VC(K)*V(I,J))*VC(K))*FY(I,J))/9.0
   ELSEIF(K.LE.8)THEN
   F(I,J,K)=F(I,J,K)*(1.0-1.0/TAO)+FEQ(I,J,K)/TAO+DT*(1.0-0.5/TAO)*((3.0*(UC(K)-U(I,J))+9.0*(UC(K)*U(I,J)+VC(K)*V(I,J))*UC(K))*FX(I,J)+(3.0*(VC(K)-V(I,J))+9.0*(UC(K)*U(I,J)+VC(K)*V(I,J))*VC(K))*FY(I,J))/36.0
   ELSE
   F(I,J,K)=F(I,J,K)*(1.0-1.0/TAO)+FEQ(I,J,K)/TAO+DT*(1.0-0.5/TAO)*((3.0*(UC(K)-U(I,J))+9.0*(UC(K)*U(I,J)+VC(K)*V(I,J))*UC(K))*FX(I,J)+(3.0*(VC(K)-V(I,J))+9.0*(UC(K)*U(I,J)+VC(K)*V(I,J))*VC(K))*FY(I,J))*4.0/9.0
   ENDIF
 ENDDO
 ENDDO
 ENDDO

END


SUBROUTINE OUTPUT
PARAMETER(IM=2401,JM=801,KM=9)
PARAMETER(ISY=400,ISB=320)
INTEGER IM,JM,KM, I,J,NSTEP
DOUBLE PRECISION X,Y
DOUBLE PRECISION U,V,D
DOUBLE PRECISION UF(IM,JM),VF(IM,JM),DF(IM,JM)
DOUBLE PRECISION DT, PI, TAO, CS
CHARACTER*8 NAME
DOUBLE PRECISION SU,SV,SX,SY,CENU,CENV,CENX,CENY,SFX,SFY,SPRKY,SXH,SYH
DOUBLE PRECISION FX,FY
DOUBLE PRECISION SUB,SVB,SXB,SYB,SFXB,SFYB,SXHB,SYHB
DOUBLE PRECISION VOR(IM,JM)
DOUBLE PRECISION ERROR, RE, DEN, U0
COMMON/IBM1/SUB(ISB),SVB(ISB),SXB(ISB),SYB(ISB),SFXB(ISB),SFYB(ISB),SXHB(ISB),SYHB(ISB)
COMMON/FORCE/FX(IM,JM),FY(IM,JM)
COMMON/IBM/SU(ISY),SV(ISY),SX(ISY),SY(ISY),CENU,CENV,CENX,CENY,SFX(ISY),SFY(ISY),SPRKY,SXH(ISY),SYH(ISY)
COMMON/STEP/NSTEP
COMMON/UVD/U(IM,JM),V(IM,JM),D(IM,JM)
COMMON/XY/X(IM,JM),Y(IM,JM)
COMMON/CONST1/DT, PI, TAO, CS
COMMON/CONST/ERROR, RE, DEN, U0

        NN=ICHAR('0')

        I0=NSTEP/100000
        I1=NSTEP/10000-I0*10
        I2=NSTEP/1000-I0*100-I1*10
        I3=NSTEP/100-I0*1000-I1*100-I2*10
        I4=NSTEP/10-I3*10-I2*100-I1*1000-I0*10000
		I5=NSTEP-I4*10-I3*100-I2*1000-I1*10000-I0*100000


        I0=NN+I0
        I1=NN+I1
        I2=NN+I2
        I3=NN+I3
        I4=NN+I4
		I5=NN+I5

        NAME='JG'//CHAR(I0)//CHAR(I1)//CHAR(I2)//CHAR(I3)//CHAR(I4)//CHAR(I5)


  DO I=1,IM
  DO J=1,JM
    UF(I,J)=U(I,J)
 	VF(I,J)=V(I,J)
 	DF(I,J)=D(I,J)
  ENDDO
  ENDDO

  DO I=2,IM-1
  DO J=2,JM-2
    VOR(I,J)=0.5*((VF(I+1,J)-VF(I-1,J))-((U(I,J+1)-U(I,J-1))))/DT
  ENDDO
  ENDDO

  DO I=2,IM-1
    VOR(I,1)=2*VOR(I,2)-VOR(I,3)
	VOR(I,JM)=2*VOR(I,JM-1)-VOR(I,JM-2)
  ENDDO

  DO J=1,JM
    VOR(1,J)=2*VOR(2,J)-VOR(3,J)
	VOR(IM,J)=2*VOR(IM-1,J)-VOR(IM-2,J)
  ENDDO

  
  OPEN(1,FILE=NAME//'.DAT')
  WRITE(1,*)'VARIABLES='
  WRITE(1,*)'"X" "Y" "D" "U" "V" "VOR"'
  WRITE(1,*) 'ZONE T="ZONE 1"'
  WRITE(1,*) 'I=',IM, ' J=',JM
  WRITE(1,*) 'F=POINT'
  DO J=1,JM
  DO I=1,IM
    WRITE(1,*) X(I,J),Y(I,J),DF(I,J),UF(I,J),VF(I,J),VOR(I,J)
  ENDDO
  ENDDO
  WRITE(1,*)'VARIABLES='
  WRITE(1,*)'"X" "Y" "D" "U" "V"'
  WRITE(1,*) 'ZONE T="ZONE 2"'
  WRITE(1,*) 'I=',ISY
  WRITE(1,*) 'F=POINT'
  DO I=1,ISY
    WRITE(1,*) SXH(I),SYH(I),DEN,U0,U0*0.0,U0*0.0
  ENDDO

  


CLOSE(1)


 ! NAME='BA'//CHAR(I0)//CHAR(I1)//CHAR(I2)//CHAR(I3)//CHAR(I4)//CHAR(I5)

  
  !OPEN(1,FILE=NAME//'.DAT')
  !WRITE(1,*)'VARIABLES='
  !WRITE(1,*)'"X" "Y"'
  !WRITE(1,*) 'ZONE T="ZONE 1"'
  !WRITE(1,*) 'I=',ISB
  !WRITE(1,*) 'F=POINT'
  
  !DO I=1,ISB
 !   WRITE(1,*) SXHB(I),SYHB(I)
 ! ENDDO
  

!CLOSE(1)



END
