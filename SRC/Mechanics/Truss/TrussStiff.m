%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     Compute the stiffness matrix of truss                       *
%*                                                                 *
%* - Call procedures:                                              *
%*     TrussStiff.m - InitTruss()                                  *
%*     TrussStiff.m - Assemble()                                   *
%*     ./ReadTruss.m - ReadTruss()                                 *
%*     SRC/Mechanics/Addres.m - Addres()                           *
%*                                                                 *
%* - Called by :                                                   *
%*     SRC/Mechanics/GetStiff.m                                    *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.21                *
%*                                                                 *
%* *****************************************************************

function TrussStiff(NG)

% Init variables of the element
InitTruss();

% Read Material and Elements
ReadTruss(NG);

fprintf('Solution phase ...\n\n');

% calculate addresses of diagonal elements
Addres(NG);

% Data check Or Solve
global cdata;
if (cdata.MODEX == 0) 
    cdata.TIM(3,:) = clock;
    cdata.TIM(4,:) = clock;
    cdata.TIM(5,:) = clock;
    return; 
end

% Assemble structure stiffness matrix
Assemble(NG);




end

% ----------------------- Functions -----------------------------------

% Init parameters of truss element
function InitTruss()
global sdata;
sdata.NNODE = 2;
sdata.NDOF = 3;

end

% Assemble structure stiffness matrix
function Assemble(NG)
global sdata;
global cdata;
S = zeros(6, 6, 'double');
ST = zeros(6, 1, 'double');
M0 = zeros(3, 1, 'double');
if( NG == 1 )
    sdata.STIFF = zeros(sdata.NWK, 1, 'double');
    if(sdata.MassType == 2)
        sdata.MASSC = zeros(sdata.NWK, 1, 'double');
    end
end


NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ; 
E = sdata.E; AREA = sdata.AREA; rho = sdata.rho; LM = sdata.LM;
for N = 1:NUME
    MTYPE = MATP(N);
    
%   compute the length of truss element
    DX = XYZ(1, N) - XYZ(4, N);
    DY = XYZ(2, N) - XYZ(5, N);
    DZ = XYZ(3, N) - XYZ(6, N);
    XL2 = DX*DX + DY*DY + DZ*DZ;
    XL = sqrt(XL2);
    
    XX = E(MTYPE) * AREA(MTYPE) * XL;
    
    ST(1) = DX / XL2;
    ST(2) = DY / XL2;
    ST(3) = DZ / XL2;
    ST(4) = -ST(1); ST(5) = -ST(2); ST(6) = -ST(3);


    for J = 1:6
        YY = ST(J) * XX;
        for I = 1:J 
            S(I, J) = ST(I)*YY; 
        end
    end
%   SRC/Mechanics/ADDBAN.m
    ADDBAN(S, LM(:, N));

    % lumped mass
    W = rho(MTYPE) * AREA(MTYPE) * XL;
    M0(1) = DX * ST(1);
    M0(2) = DY * ST(2);
    M0(3) = DZ * ST(3);

    a = 1/2/sum(M0);
    M = a * W * [M0;M0];
    ADDMASS(M, LM(:, N));

    % Compatible Mass Matrix
    if(sdata.MassType == 2)
        R = [DX DY DZ 0 0 0;0 0 0 DX DY DZ]/XL;
        Me = [1/3 1/6;1/6 1/3];
        M0 = W*R'*Me*R;
        ADDMASSC(M0, LM(:, N));
    end

    
end

% The third time stamp
cdata.TIM(3, :) = clock;

end
