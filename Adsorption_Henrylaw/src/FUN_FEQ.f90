FUNCTION FUN_FEQ(K1,u1,v1,C_INF)
    USE SIMPLE
    IMPLICIT NONE
    
    INTEGER K1
    REAL FUN_FEQ,u1,V1,C_INF
    
    FUN_FEQ=W(K1)*c_INF*(1.0+3.0*(CX(K1)*U1+CY(K)*V1/CK)+4.5*(U1*cx(k)+V1*cy(k))**2-1.5*(U1**2+V1**2))
    
RETURN
END