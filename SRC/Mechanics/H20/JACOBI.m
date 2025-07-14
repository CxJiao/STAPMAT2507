%* *****************************************************************
%* - Function of STAPMAT in stiffness and solver phase             *
%*                                                                 *
%* - Purpose:                                                      *
%*     To calculate jacobi matrix at point(xi,eta,zeta)            *
%*                                                                 *
%* - Call procedures: None                                         *
%*                                                                 *
%* - Called by :                                                   *
%*     ./ReadH20.m - ReadAssembleLoads()                           *
%*     ./H20Stiff.m - Assemble()                                   *
%*     ./B_cal()                                                   *
%*                                                                 *
%* *****************************************************************

function jaco = JACOBI(xi,eta,zeta,N)

global sdata;
XYZ = sdata.XYZ2(:,N);
xyz = (reshape(XYZ,3,20))';

jaco = dN_cal(xi,eta,zeta) * xyz;

end