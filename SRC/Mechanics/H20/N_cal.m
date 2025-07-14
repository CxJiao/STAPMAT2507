%* *****************************************************************
%* - Function of STAPMAT in stiffness and solver phase             *
%*                                                                 *
%* - Purpose:                                                      *
%*     To calculate shape function                                 *
%*                                                                 *
%* - Call procedures: None                                         *
%*                                                                 *
%* - Called by :                                                   *
%*     ./ReadH20.m - ReadAssembleLoads()                           *
%*     ./H20Stiff.m - Assemble()                                   *
%*     ./H20Stressf.m                                              *
%*                                                                 *
%* *****************************************************************


function NN = N_cal(xi,eta,zeta)
% Node coordinate(LOCAL)
XI = [-1 1 1 -1 -1 1 1 -1 0 0 0 0 1 1 -1 -1 -1 1 1 -1];
ETA = [-1 -1 1 1 -1 -1 1 1 -1 1 1 -1 0 0 0 0 -1 -1 1 1];
ZETA = [-1 -1 -1 -1 1 1 1 1 -1 -1 1 1 -1 1 1 -1 0 0 0 0];

NN = zeros(3,60);
N = zeros(20,1);% diagonal elements

for i = 1:8
    xi0 = XI(i)*xi; eta0 = ETA(i)*eta; zeta0 = ZETA(i)*zeta;
    Ni = (1+xi0) * (1+eta0) * (1+zeta0) * (xi0+eta0+zeta0-2)/8;
    N(i) = Ni;
end

for i = 9:12
    eta0 = ETA(i)*eta; zeta0 = ZETA(i)*zeta;
    Ni = (1-xi^2) * (1+eta0) * (1+zeta0)/4;
    N(i) = Ni;
end

for i = 13:16
    xi0 = XI(i)*xi; zeta0 = ZETA(i)*zeta;
    Ni = (1-eta^2) * (1+xi0) * (1+zeta0)/4;
    N(i) = Ni;
end

for i = 17:20
    xi0 = XI(i)*xi; eta0 = ETA(i)*eta; 
    Ni = (1-zeta^2) * (1+xi0) * (1+eta0)/4;
    N(i) = Ni;
end

for i = 1:20
    NN(1,3*i-2)=N(i); NN(2,3*i-1)=N(i); NN(3,3*i)=N(i);
end

end