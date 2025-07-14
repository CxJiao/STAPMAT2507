%* *****************************************************************
%* - Function of STAPMAT in stiffness and solver phase             *
%*                                                                 *
%* - Purpose:                                                      *
%*     To calculate geometry matrix in 3D                          *
%*                                                                 *
%* - Call procedures: None                                         *
%*                                                                 *
%* - Called by :                                                   *
%*     ./H20Stiff.m - Assemble()                                   *
%*     ./H20Stressf.m                                              *
%*                                                                 *
%* *****************************************************************


function B = B_cal(xi,eta,zeta,N)

jaco = JACOBI(xi,eta,zeta,N);
invJ = inv(jaco);
dN = invJ * dN_cal(xi,eta,zeta);


B = zeros(6,60);
for i =1:20
    B(1,3*i-2) = dN(1,i);
    B(2,3*i-1) = dN(2,i);
    B(3,3*i) = dN(3,i);
    B(4,3*i-2) = dN(2,i); B(4,3*i-1) = dN(1,i);
    B(5,3*i-1) = dN(3,i); B(5,3*i) = dN(2,i);
    B(6,3*i-2) = dN(3,i); B(6,3*i) = dN(1,i);
end

end