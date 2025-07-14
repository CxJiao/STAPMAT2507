%* *****************************************************************
%* - Function of STAPMAT in stiffness and solver phase             *
%*                                                                 *
%* - Purpose:                                                      *
%*     To calculate Derivative of shape function                   *
%*                                                                 *
%* - Call procedures: None                                         *
%*                                                                 *
%* - Called by :                                                   *
%*     ./Jacobi()                                                  *
%*     ./B_cal()                                                   *
%*                                                                 *
%* *****************************************************************


function dN = dN_cal(xi,eta,zeta)
% Node coordinate(LOCAL)
XI = [-1 1 1 -1 -1 1 1 -1 0 0 0 0 1 1 -1 -1 -1 1 1 -1];
ETA = [-1 -1 1 1 -1 -1 1 1 -1 1 1 -1 0 0 0 0 -1 -1 1 1];
ZETA = [-1 -1 -1 -1 1 1 1 1 -1 -1 1 1 -1 1 1 -1 0 0 0 0];

dN = zeros(3,20);

for i = 1:8
    xi0 = XI(i)*xi; eta0 = ETA(i)*eta; zeta0 = ZETA(i)*zeta;
    dN(1,i) = (1+eta0) * (1+zeta0) * (2*xi0+eta0+zeta0-1) * XI(i)/8;
    dN(2,i) = (1+xi0) * (1+zeta0) * (xi0+2*eta0+zeta0-1) * ETA(i)/8;
    dN(3,i) = (1+xi0) * (1+eta0) * (xi0+eta0+2*zeta0-1) * ZETA(i)/8;
end

for i = 9:12
    eta0 = ETA(i)*eta; zeta0 = ZETA(i)*zeta;
    dN(1,i) = -2 * xi * (1+eta0) * (1+zeta0)/4;
    dN(2,i) = (1-xi^2) * ETA(i) * (1+zeta0)/4;
    dN(3,i) = (1-xi^2) * (1+eta0) * ZETA(i)/4;
end

for i = 13:16
    xi0 = XI(i)*xi; zeta0 = ZETA(i)*zeta;
    dN(1,i) = XI(i) * (1-eta^2) * (1+zeta0)/4;
    dN(2,i) = -2 * (1+xi0) * eta * (1+zeta0)/4;
    dN(3,i) = (1+xi0) * (1-eta^2) * ZETA(i)/4;
end

for i = 17:20
    xi0 = XI(i)*xi; eta0 = ETA(i)*eta;
    dN(1,i) = XI(i) * (1+eta0) * (1-zeta^2)/4;
    dN(2,i) = (1+xi0) * ETA(i) * (1-zeta^2)/4;
    dN(3,i) = -2 * (1+xi0) * (1+eta0) * zeta/4;
end

end