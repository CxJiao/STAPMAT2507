%* *****************************************************************
%* - Function of STAPMAT in Solver phase                           *
%*                                                                 *
%* - Purpose:                                                      *
%*     To solve finite element dynamic equilibrium equations       *
%*                                                                 *
%* - Call procedures:                                              *
%*     Dynamic.m                 - KMul()                          *
%*     Dynamic.m                 - WriteDVA()                      *
%*     Dynamic.m                 - DLDLT()                         *
%*     Dynamic.m                 - EquSol()                        *
%*                                                                 *
%* - Called by :                                                   *
%*     solve.m                                                     *
%*                                                                 *
%* *****************************************************************
function Dynamics(NLC)
% input: NLC - number of load case

global sdata;
ID = sdata.ID;
InitialDis = sdata.InitialDis;
pM = sdata.CparaM; pK = sdata.CparaK;
Ddis = sdata.Ddis; Dvol = sdata.Dvol; Dacc = sdata.Dacc; 
af = sdata.DyAlphaF; am = sdata.DyAlphaM; g = sdata.DyGamma; b = sdata.DyBeta;
ENDT = sdata.ENDTIME; dt = sdata.TimeStep; NSTEPS = sdata.NSTEPS;
sdata.EnKine = zeros(NSTEPS,1); EnKine = sdata.EnKine;
NVNL = sdata.NVNL; DLDC = sdata.DLDC; DNOD = sdata.DNOD;
DDIRE = sdata.DDIRE; DVNL = sdata.DVNL; R = sdata.R;

STIFF = sdata.STIFF; MAXA = sdata.MAXA;
if(sdata.MassType == 1)
    MASS = sdata.MASS;
elseif(sdata.MassType == 2)
    MASS = sdata.MASSC;
end

global cdata;
IOUT = cdata.IOUT;
fprintf(IOUT, '\n\n D Y N A M I C  S O L U T I O N S: LOAD CASE %3d', NLC);

% Set Parameter Default Value
if(af == 0 && am == 0 && g == 0 && b == 0)
    af = 1/2; am = 1/2; b = 1/4; g = 1/2;
    sdata.DyAlphaF = af; sdata.DyAlphaM = am; sdata.DyGamma = g; sdata.DyBeta = b;
end

% check unconditional stability
if (~(af>=0 && af<=1 && am>=0 && am<=1 && g>=0 && g<=3/2 && b>=0 && b<=1))
    error(' *** ERROR *** DYNAMICS PARAMETERS ARE OUT OF RANGE');
end
if(~(am<=af && af<=1/2 && b>=1/4+(af-am)/2))
    warning(' *** WARNING *** Algorithm May Be Conditional Stable');
end

% Static Solution AS Initial Displacements
if(InitialDis == 2)
    STIFF_LDLT = DLDLT(STIFF);
    Ddis = EquSol(STIFF_LDLT,R(:,NLC));
    sdata.Ddis = Ddis;
end

WriteDVA(0);
if(sdata.MassType == 1)
    EnKine(1) = sum(MASS.*(Dvol.^2))/2;
elseif(sdata.MassType == 2)
    EnKine(1) = Dvol'*KMul(MASS,Dvol)/2;
end


% Generalized alpha method Parameters
bdt = b*dt; bdt2 = bdt*dt;
ck = 1-af; c0 = (1-am)/bdt2; c1 = ck*g/bdt;
c2 = dt*c0; c3 = c2*dt/2-1; c4 = ck*g/b-1; c5 = ck*(g/b/2-1)*dt;
c1a = c1*pM; c4a = c4*pM; c5a = c5*pM;
c1b = c1*pK; c4b = c4*pK; c5b = c5*pK;

K_hat = (ck+c1b)*STIFF; M_hat = (c0+c1a)*MASS;
if(sdata.MassType == 1)
    for I = 1:length(MASS)
        K_hat(MAXA(I)) = K_hat(MAXA(I))+M_hat(I);
    end
elseif(sdata.MassType == 2)
    K_hat = K_hat+M_hat;
end
K_hatF = DLDLT(K_hat);

Sf = R(:,NLC); Df1 = zeros(size(Sf)); Df2 = zeros(size(Sf));
% Initialize Time-Variant Load
NLCloc = [];% Time-Variant Load in this load case
for I = 1:NVNL
    if(DLDC(I) == NLC)
        NLCloc = [NLCloc I];
        INDEX = ID(DDIRE(I),DNOD(I));
        Df2(INDEX) = DVNL(1,I);
    end
end


for I = 2:NSTEPS
    Df1 = Df2; Df2 = zeros(size(Sf));
    Ddis = sdata.Ddis; Dvol = sdata.Dvol; Dacc = sdata.Dacc;
    for J = 1:length(NLCloc)
            INDEX = ID(DDIRE(NLCloc(J)),DNOD(NLCloc(J)));
            Df2(INDEX) = DVNL(I,NLCloc(J));
    end
    Q = Sf + (1-af)*Df2 + af*Df1;
    
    Kfac = (c1b-af)*Ddis + c4b*Dvol + c5b*Dacc;
    Mfac = (c0+c1a)*Ddis + (c2+c4a)*Dvol + (c3+c5a)*Dacc;
    
    if(sdata.MassType == 1)
        Q_hat = Q + KMul(STIFF,Kfac) + MASS.*Mfac;
    elseif(sdata.MassType == 2)
        Q_hat = Q + KMul(STIFF,Kfac) + KMul(MASS,Mfac);
    end

    % solve the equation
    Ddis2 = EquSol(K_hatF,Q_hat);
    Dvol2 = (Ddis2-Ddis)*g/bdt+(1-g/b)*Dvol+(1-g/2/b)*dt*Dacc;
    Dacc2 = (Ddis2-Ddis)/bdt2-Dvol/bdt-(1/2/b-1)*Dacc;

    sdata.Ddis = Ddis2; sdata.Dvol = Dvol2; sdata.Dacc = Dacc2;
    WriteDVA(min((I-1)*dt,ENDT));% print current results
    if(sdata.MassType == 1)
        EnKine(I) = sum(MASS.*(Dvol2.^2))/2;
    elseif(sdata.MassType == 2)
        EnKine(I) = Dvol2'*KMul(MASS,Dvol2)/2;
    end
end
% print Kinetic Energy
TimeStamp = 0:dt:ENDT;
if(TimeStamp(end) ~= ENDT) TimeStamp = [TimeStamp ENDT]; end
fprintf(IOUT, '\n\n KINETIC ENERGY FOR LOAD CASE %3d', NLC);
fprintf(IOUT, '\n TIME               '); fprintf(IOUT,'%18.6e',TimeStamp);
fprintf(IOUT, '\n KINETIC_ENERGT     '); fprintf(IOUT,'%18.6e',EnKine);

end


% Matrix Multiplication for One-Dimensional Bandwidth-Variant STIFFNESS
function result = KMul(STIFF,x)

global sdata;
MHT = sdata.MHT; MAXA = sdata.MAXA;
% STIFF = sdata.STIFF;

result = zeros(size(x));
for ind = 1:length(MHT)
    % column
    ed = ind; st = ind-MHT(ind);
    result(st:ed) = result(st:ed) + flip(STIFF(MAXA(ind):MAXA(ind+1)-1))*x(ind);

    % row
    if(st<ed)
        result(ind) = result(ind) + sum(flip(STIFF(MAXA(ind)+1:MAXA(ind+1)-1)).*x(st:ed-1));
    end

end
end


% Print Displacements & Volecity & Acceleration
function WriteDVA(Time)

% Get global data
global cdata;
global sdata;
IOUT = cdata.IOUT;
NUMNP = cdata.NUMNP;
ID = sdata.ID;
Ddis = sdata.Ddis; Dvol = sdata.Dvol; Dacc = sdata.Dacc;

fprintf(IOUT, '\n\n Time = %3d', Time);


fprintf(IOUT, ['\n\n D I S P L A C E M E N T S\n' ...
    '\n         NODE           X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT\n']);

D = zeros(3, 1, 'double');
for II = 1:NUMNP
    D(:) = 0;
    if (ID(1, II) ~= 0) D(1) = Ddis(ID(1, II)); end
    if (ID(2, II) ~= 0) D(2) = Ddis(ID(2, II)); end
    if (ID(3, II) ~= 0) D(3) = Ddis(ID(3, II)); end

    fprintf(IOUT, ' %10d        %18.6e%18.6e%18.6e\n', II, D(1), D(2), D(3));
end


fprintf(IOUT, ['\n\n V O L E C I T I E S\n' ...
    '\n         NODE           X-VOLECITY        Y-VOLECITY        Z-VOLECITY\n']);
V = zeros(3, 1, 'double');
for II = 1:NUMNP
    V(:) = 0;
    if (ID(1, II) ~= 0) V(1) = Dvol(ID(1, II)); end
    if (ID(2, II) ~= 0) V(2) = Dvol(ID(2, II)); end
    if (ID(3, II) ~= 0) V(3) = Dvol(ID(3, II)); end

    fprintf(IOUT, ' %10d        %18.6e%18.6e%18.6e\n', II, V(1), V(2), V(3));
end


fprintf(IOUT, ['\n\n A C C E L E R A T I O N S\n' ...
    '\n         NODE           X-ACCELERATION    Y-ACCELERATION    Z-ACCELERATION\n']);
A = zeros(3, 1, 'double');
for II = 1:NUMNP
    A(:) = 0;
    if (ID(1, II) ~= 0) A(1) = Dacc(ID(1, II)); end
    if (ID(2, II) ~= 0) A(2) = Dacc(ID(2, II)); end
    if (ID(3, II) ~= 0) A(3) = Dacc(ID(3, II)); end

    fprintf(IOUT, ' %10d        %18.6e%18.6e%18.6e\n', II, A(1), A(2), A(3));
end
end

% LDLT factorization for K_hat
function A = DLDLT(A)

% Get global data
global sdata;
MAXA = sdata.MAXA; NEQ = sdata.NEQ; 

for N = 1:NEQ
    KN = MAXA(N);
    KL = KN + 1;
    KU = MAXA(N + 1) - 1;
    KH = KU - KL;
    
    if (KH > 0)
        K = N - KH;
        IC = 0;
        KLT = KU;
        for J = 1:KH
            IC = IC + 1;
            KLT = KLT - 1;
            KI = MAXA(K);
            ND = MAXA(K+1) - KI - 1;
            if (ND > 0)
                KK = min(IC, ND);
                C = 0.0;
                for L = 1:KK C = C+A(KI+L)*A(KLT+L); end
                A(KLT) = A(KLT) - C;
            end
            K = K + 1;
        end
    end
    
    if (KH >= 0)
        K = N;
        B = 0.0;
        for KK = KL:KU
            K = K - 1;
            KI = MAXA(K);
            C = A(KK) / A(KI);
            B = B + C*A(KK);
            A(KK) = C;
        end
        A(KN) = A(KN) - B;
    end
    
    if (A(KN) <= 0)
        error(['STOP - (Dynamic Solution) Stiffness matrix is not positive definite\n' ...
            'Nonpositive number for equation %8d is %20.12e\n'], N, A(KN));
    end
end

end


% solve equition
function result = EquSol(A,R)

% Get global data
global sdata;
MAXA = sdata.MAXA; NEQ = sdata.NEQ;

% Reduce right-hand-side load vector
for N = 1:NEQ
    KL = MAXA(N) + 1;
    KU = MAXA(N+1) - 1;
    if (KU-KL >= 0)
        K = N;
        C = 0.0;
        for KK = KL:KU
            K = K - 1;
            C = C + A(KK) * R(K);
        end
        R(N) = R(N) - C;
    end
end

% Back-Substitute
for N = 1:NEQ
    K = MAXA(N);
    R(N) = R(N) / A(K);
end

if (NEQ == 1) return; end;

N = NEQ;
for L = 2:NEQ
    KL = MAXA(N) + 1;
    KU = MAXA(N+1) - 1;
    if (KU-KL >= 0)
        K = N;
        for KK = KL:KU
            K = K - 1;
            R(K) = R(K) - A(KK)*R(N);
        end
    end
    N = N - 1;
end

result = R(:);

end