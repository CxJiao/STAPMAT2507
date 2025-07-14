%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     Compute the stiffness matrix of H20                         *
%*                                                                 *
%* - Call procedures:                                              *
%*     H20Stiff.m - InitH20()                                      *
%*     H20Stiff.m - Assembel()                                     *
%*     ./ReadH20.m - ReadH20()                                     *
%*     SRC/Mechanics/Addres.m - Addres()                           *
%*                                                                 *
%* - Called by :                                                   *
%*     SRC/Mechanics/GetStiff.m                                    *
%*                                                                 *
%* *****************************************************************

function H20Stiff(NG)

% Init variables of the element
InitH20();

% Read Material and Elements
ReadH20(NG);

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
function InitH20()
global sdata;
sdata.NNODE = 20;
sdata.NDOF = 3;

end



% Assemble structure stiffness matrix
function Assemble(NG)
global sdata;
global cdata;

if( NG == 1 )
    sdata.STIFF = zeros(sdata.NWK, 1, 'double');
    if(sdata.MassType == 2)
        sdata.MASSC = zeros(sdata.NWK, 1, 'double');
    end
end

NUME = sdata.NUME2; MATP = sdata.MATP2;
E = sdata.E2;nu = sdata.nu; rho = sdata.rho; LM = sdata.LM2;
GC3 = sdata.GC3; GW3 = sdata.GW3;

for N = 1:NUME
    MTYPE = MATP(N); v = nu(MTYPE);

    % Physics Matrix
    D = E(MTYPE) / (1 + v) / (1 - 2*v) * [1-v v v 0 0 0;...
        v 1-v v 0 0 0;v v 1-v 0 0 0;0 0 0 (1-2*v)/2 0 0;...
        0 0 0 0 (1-2*v)/2 0;0 0 0 0 0 (1-2*v)/2];

    % Element Stiffness Matrix
    S=zeros(60,60);

    for i = 1:3
        for j = 1:3
            for k = 1:3
            S = S + GW3(i)*GW3(j)*GW3(k) * B_cal(GC3(i),GC3(j),GC3(k),N)' * D *...
            B_cal(GC3(i),GC3(j),GC3(k),N) * det(JACOBI(GC3(i),GC3(j),GC3(k),N));
            end
        end
    end
%   SRC/Mechanics/ADDBAN.m
    ADDBAN(S, LM(:, N));

    
    Volume = 0;
    for i = 1:3
        for j = 1:3
            for k = 1:3
                Volume = Volume + GW3(i)*GW3(j)*GW3(k) * det(JACOBI(GC3(i),GC3(j),GC3(k),N));
            end
        end
    end
    W = rho(MTYPE) * Volume;
    
    M0 = zeros(60,60);
    for i = 1:3
        for j = 1:3
            for k = 1:3
                M0 = M0 + GW3(i)*GW3(j)*GW3(k) * N_cal(GC3(i),GC3(j),GC3(k))'...
                    * N_cal(GC3(i),GC3(j),GC3(k)) * det(JACOBI(GC3(i),GC3(j),GC3(k),N));
            end
        end
    end

    % Lumped Mass Matrix
    a = 1/sum(diag(M0));
    M = (3*a) * W * diag(M0);
    ADDMASS(M, LM(:, N));

    % Compatible Mass Matrix
    if(sdata.MassType == 2)
        M0 = rho(MTYPE) * M0;
        ADDMASSC(M0, LM(:, N));
    end

end

% The third time stamp
cdata.TIM(3, :) = clock;

end