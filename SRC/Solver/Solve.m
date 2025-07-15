%* *****************************************************************
%* - Function of STAPMAT in Solver phase                           *
%*                                                                 *
%* - Purpose:                                                      *
%*     To solve finite element static equilibrium equations        *
%*                                                                 *
%* - Call procedures:                                              *
%*     ./LDLTFactor.m            - LDLTFactor()                    *
%*     Solve.m                   - Stiff2Sparse()                  *
%*     ./ColSol.m                - ColSol()                        *  
%*     Solve.m                   - WriteDis()                      *
%*     SRC/Mechanics/GetStress.m - GetStress()                     *
%*                                                                 *
%* - Called by :                                                   *
%*     stapmat.m                                                   *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.22                *
%*                                                                 *
%* *****************************************************************

function Solve()

global cdata;
global sdata;

NEQ = sdata.NEQ;
NLCASE = cdata.NLCASE;
MODEX = cdata.MODEX;
sdata.DIS = zeros(NEQ, NLCASE, 'double');
sdata.STRAIN = zeros(NEQ, NLCASE, 'double');
sdata.STRESS = zeros(NEQ, NLCASE, 'double');

% The pre-process of Solution
if (MODEX == 2)%%2或3时启动动力学求解程序，其中2为隐式的广义α法

    for L = 1:NLCASE Dynamics(L); end
    cdata.TIM(4,:) = clock;
    
elseif (MODEX == 3)%%3为显式的速度verlet法
    for L = 1:NLCASE Dynamicsverlet(L); end
    cdata.TIM(4,:) = clock;
    
else

    if (MODEX == 1) LDLTFactor();
    else SPSTIFF = Stiff2Sparse(); end

    cdata.TIM(4,:) = clock;

    % Solve
    for L = 1:NLCASE

        %   Solve the equilibrium equations to calculate the displacements
        if (MODEX == 1) ColSol(L);
        else sdata.DIS(:,L) = SPSTIFF \ sdata.R(:,L); end

        %   Print displacements
        WriteDis(L);

        %   Calculation of stresses
        GetStress(L);

    end
end

cdata.TIM(5, :) = clock;

end

% ----------------------- Functions -----------------------------------

% Convert the stiff vector to a sparse stiff matrix
function SPSTIFF = Stiff2Sparse()

global sdata;
A = sdata.STIFF; MAXA = sdata.MAXA; NEQ = sdata.NEQ; NWK = sdata.NWK;

IIndex = zeros(NWK*2-NEQ, 1);
JIndex = IIndex;
STIFF = IIndex;
NUM = 1;
NUMC = 0;
for N = 1:NEQ
    KU = MAXA(N + 1) - MAXA(N);
    for L = 1:KU
        IIndex(NUM) = N;
        JIndex(NUM) = N - L + 1;
        STIFF(NUM) = A(NUM);
        NUM = NUM + 1;
        if (L == 1) NUMC = NUMC + 1;continue; end
        SYMN = NUM-1 - NUMC + NWK;
        IIndex(SYMN) = N - L + 1;
        JIndex(SYMN) = N;
        STIFF(SYMN) = A(NUM-1);
    end
end

SPSTIFF = sparse(IIndex, JIndex, STIFF, NEQ, NEQ);
end

% Print Displacements
function WriteDis(NUM)

% Get global data
global cdata;
global sdata;
IOUT = cdata.IOUT;
NUMNP = cdata.NUMNP;
DIS = sdata.DIS(:, NUM); ID = sdata.ID;

fprintf(IOUT, '\n\n LOAD CASE %3d', NUM);
fprintf(IOUT, ['\n\n D I S P L A C E M E N T S\n' ...
    '\n       NODE           X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT\n']);

D = zeros(3, 1, 'double');
for II = 1:NUMNP
    D(:) = 0;
    if (ID(1, II) ~= 0) D(1) = DIS(ID(1, II)); end
    if (ID(2, II) ~= 0) D(2) = DIS(ID(2, II)); end
    if (ID(3, II) ~= 0) D(3) = DIS(ID(3, II)); end
    
    fprintf(IOUT, ' %10d        %18.6e%18.6e%18.6e\n', II, D(1), D(2), D(3));
end

end