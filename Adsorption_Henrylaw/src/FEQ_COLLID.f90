SUBROUTINE FEQ_COLLID
    USE SIMPLE
    IMPLICIT NONE
    real q_k(0:8,0:NX,0:NY)
    
    do j=0,NY
    do i=0,NX
	do k=0,8
        feq(k)=w(k)*C(i,j)*(1.0+3.0*(u(i,j)*cx(k)+v(i,j)*cy(k))+4.5*(u(i,j)*cx(k)+v(i,j)*cy(k))**2-1.5*(U(I,J)**2+V(I,J)**2))
        q_k(k,i,j)=W(K)*(1-0.5*OMEGA)*QXY(I,J)!*(1+3.0*(u(i,j)*cx(k)+v(i,j)*cy(k)))
		f(k,i,j)=omega*feq(k)+(1.0-omega)*f(k,i,j)+DT*q_k(k,i,j)
	end do
	end do
    end do
RETURN
END