﻿!=======simulation practice 1st========!
!====channel flow====!
!inlet: u=; outlet:fully developed; upper&lower: non-slip
PROGRAM MAIN
implicit none
integer,parameter::NX=4001,NY=101		!size of simulation domain
integer,parameter::NS=200				!number of marker nodes on the boundary
integer,parameter::QD=9
integer,parameter::T_plot=300000
integer,parameter::T_MAX =300000
double precision f(NX,NY,QD),feq(NX,NY,QD),rho(NX,NY),U(NX,NY),V(NX,NY),Fix(NX,NY),Fiy(NX,NY),ff(NX,NY,QD)
double precision x(NX,NY),y(NX,NY),U0,Re,DEN,gradP
double precision cx(QD),cy(QD),w(QD),tau,omega,niu,c,cs,opp(QD)
double precision u_analy(NY),ERROR
real start, finish

double precision lx,ly,dx,dy,dt,pi,temp1,temp2,DD,UU,VV
integer i,j,k,i1,i2,j1,j2,M,N,NSTEP,I0,I3,I4,I5,ZZ
CHARACTER*8 NAME,NAME1 

!double precision xs(NS),ys(NS),cenx,ceny,ds,Rs,D
!double precision Us(NS),Vs(NS),Fsx(NS),Fsy(NS)
!double precision CD,CL,belta
lx=40.0
ly=1.0
dx=lx/(NX-1)
dy=ly/(NY-1)
dt=dx
c=dx/dt
cs=c/sqrt(3.0)
pi=4.0*atan(1.0)
U0=0.1
Re=100.0
DEN=1.0
niu=U0*ly/Re
tau=niu/(cs**2*dt)+0.5
omega=1.0/tau
opp(1:9)=(/3,4,1,2,7,8,5,6,9/)
do i=1,4
    cx(i)=cos(0.5*pi*(i-1))*dx/dt
    cy(i)=sin(0.5*pi*(i-1))*dx/dt
    w(i)=1.0/9.0
enddo
do i=5,8
    cx(i)=cos(0.5*pi*(i-9.0/2.0))*sqrt(2.0)*dx/dt
    cy(i)=sin(0.5*pi*(i-9.0/2.0))*sqrt(2.0)*dx/dt
    w(i)=1.0/36.0
enddo
    cx(9)=0.0
    cy(9)=0.0
    w(9)=4.0/9.0
!coordinates of the Eularian points
do i=1,NX
do j=1,NY
x(i,j)=lx*(i-1)/(nx-1)
y(i,j)=ly*(j-1)/(ny-1)
enddo
enddo
U   =   U0
V	=   0.0
rho =   DEN
do i=1,nx
do j=1,ny
	do k=1,QD
        temp1=cx(k)*U(i,j)+cy(k)*V(i,j)
        temp2=U(i,j)**2+V(i,j)**2
        feq(i,j,k)=w(k)*(rho(i,j)+3.0*temp1+4.5*temp1*temp1-1.5*temp2) 
	enddo
enddo
enddo
f=feq
!==========================================!
!ANALYTICAL SOLUTION OF THE fully developed
!VELOCITY PROFILE
!Based on the power-law fluid flow in a channel
!Ref.(F.B. Tian,2014CM; R. P. Chhabra, 2008)
do j=1,ny
u_analy(j)=1.5*U0*(1.0-(1-2.0*y(1,j)/ly)**2)
enddo
ERROR=1.0
NSTEP=0
OPEN(11,FILE='LOG.PLT')
WRITE(11,*) 'variables=Timestep, U/U0, ERROR'
!!Loop begin
Call CPU_TIME(START)
print*, start
DO WHILE((NSTEP.LE.T_MAX).OR.(ERROR.LT.1.0D-3))
!!Rasidual Error    
DO J=1,NY
    TEMP1=TEMP1+ABS(U(NX/2,J)-U_ANALY(J))
ENDDO
ERROR=TEMP1/NY
    
!!The following are LBM processes
WRITE(11,*)NSTEP,U((NX-1)/2,(NY-1)/2)/U0,ERROR
WRITE(*,*)NSTEP,U((NX-1)/2,(NY-1)/2)/U0,ERROR
 NSTEP=NSTEP+1
 
!========FEQ=collid======!
do i=1,nx
do j=1,ny
do k=1,QD
    temp1=cx(k)*U(i,j)+cy(k)*V(i,j)
    temp2=U(i,j)**2+V(i,j)**2
    feq(i,j,k)=w(k)*rho(i,j)*(1.0+3.0*temp1+4.5*temp1*temp1-1.5*temp2)
	f(i,j,k)=(1-omega)*f(i,j,k)+omega*feq(i,j,k) 
enddo
enddo
enddo
!======stream=======!
ff=f
do i=2,nx
do j=1,ny
	f(i,j,1)=ff(i-1,j,1)
enddo
enddo
do i=1,nx
do j=2,ny
	f(i,j,2)=ff(i,j-1,2)
enddo
enddo
do i=1,nx-1
do j=1,ny
	f(i,j,3)=ff(i+1,j,3)
enddo
enddo
do i=1,nx
do j=1,ny-1
	f(i,j,4)=ff(i,j+1,4)
enddo
enddo
do i=2,nx
do j=2,ny
	f(i,j,5)=ff(i-1,j-1,5)
enddo
enddo
do i=1,nx-1
do j=2,ny
	f(i,j,6)=ff(i+1,j-1,6)
enddo
enddo
do i=1,nx-1
do j=1,ny-1
	f(i,j,7)=ff(i+1,j+1,7)
enddo
enddo
do i=2,nx
do j=1,ny-1
	f(i,j,8)=ff(i-1,j+1,8)
enddo
enddo
!===Macro-parameter===!
do i=1,nx
do j=1,ny  
dd=0.0
uu=0.0
vv=0.0
do k=1,QD
	dd=dd+f(i,j,k)
	uu=uu+f(i,j,k)*cx(k)
	vv=vv+f(i,j,k)*cy(k)
enddo
rho(i,j)=   dd
u(i,j)= uu/dd !+0.5*Fix(i,j)*dt
v(i,j)= vv/dd !+0.5*Fiy(i,j)*dt
enddo
enddo
DO J=1,NY
    U(1,J)   =  U0
    V(1,J)   =  0.0
    RHO(1,J) =  RHO(2,J)
ENDDO
DO I=1,NX
    U(I,NY) =0.0
    V(I,NY) =0.0
    U(I,1)  =0.0
    V(I,1)  =0.0
    RHO(I,NY)=RHO(I,NY-1)
    RHO(I,1) =RHO(I,2)
ENDDO
DO J=1,NY
    U(NX,J)  = U0!U(NX-1,J)
    V(NX,J)  = 0.0
    RHO(NX,J)= RHO(NX-1,J)
ENDDO


!==boundary condition==
!non-equilibrium extrapolation
DO K=1,QD
DO I=1,NX
    F(I,1,k)    =FEQ(I,1,k) +F(I,2,opp(k))   -FEQ(I,2,opp(k))
    F(I,NY,k)   =FEQ(I,NY,k)+F(I,NY-1,opp(k))-FEQ(I,NY-1,opp(k))
ENDDO
DO J=1,NY
    F(1,J,k)    =FEQ(1,J,k) +F(2,J,opp(k))   -FEQ(2,J,opp(k))
    F(NX,J,k)   =FEQ(NX,J,k)+F(NX-1,J,opp(k))-FEQ(NX-1,J,opp(k))
ENDDO
ENDDO 
! density at corner nodes - 2nd extrapolation from bulk nodes
rho(1,1)     = 2*rho(2,2)        - rho(3,3)
rho(NX,1)    = 2*rho(NX-1,2)     - rho(NX-2,3)
rho(NX,NY)   = 2*rho(NX-1,NY-1)  - rho(NX-2,NY-2)
rho(1,NY)    = 2*rho(2,NY-1)     - rho(3,NY-2)

!!==========================================================!
!!after t_plot time output the field data to tecplot files
!!---------------------------------------------------------
	if(mod(NStep,T_plot).EQ.0)then
		ZZ=ICHAR('0')
        I0=NSTEP/100000
	    I1=NSTEP/10000-I0*10
	    I2=NSTEP/1000-I0*100-I1*10
	    I3=NSTEP/100-I0*1000-I1*100-I2*10
	    I4=NSTEP/10-I3*10-I2*100-I1*1000-I0*10000
	    I5=NSTEP-I4*10-I3*100-I2*1000-I1*10000-I0*100000
	    I0=ZZ+I0
	    I1=ZZ+I1
	    I2=ZZ+I2
	    I3=ZZ+I3
	    I4=ZZ+I4
	    I5=ZZ+I5
        NAME=CHAR(I0)//CHAR(I1)//CHAR(I2)//CHAR(I3)//CHAR(I4)//CHAR(I5)
        
		OPEN(1,FILE='OUTPUT//F'//NAME//'.PLT')
		WRITE(1,*)'VARIABLES='
		WRITE(1,*)'"X" "Y" "U" "V" "D" '
		WRITE(1,*) 'ZONE T="ZONE 1"'
		WRITE(1,*) 'I=',NX, ' J=',NY
		WRITE(1,*) 'F=POINT'
		DO J=1,NY
			DO I=1,NX
			WRITE(1,*) X(I,J),Y(I,J),U(I,J),V(I,J),rho(I,J)
			ENDDO
		ENDDO
		CLOSE(1)
	ELSE
	ENDIF

enddo

    OPEN(2,FILE='OUTPUT//RAT_P'//NAME//'.PLT')
    do J=1,NY
        TEMP1=U((NX-1)/2,J)/U0
    write(2,*)Y((NX-1)/2,J),TEMP1
    enddo
    close(2)
    OPEN(3,FILE='OUTPUT//RE_P'//NAME//'.PLT')
    do J=1,NY
    write(3,*)Y((NX-1)/2,J),U((NX-1)/2,J)
    enddo
    close(3)
    
!==================================================！   
    Call CPU_TIME(FINISH)
    PRINT*,'The program start at:', START
    PRINT*,'The program FINISH at:', FINISH
    
    
    OPEN(1,FILE='OUTPUT//OUT.TXT')
    WRITE(1,*) 'THE PROGRAM COST TOTAL TIME:'
    WRITE(1,*) FINISH-START
END PROGRAM