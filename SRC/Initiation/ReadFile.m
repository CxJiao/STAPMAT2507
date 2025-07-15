%* *****************************************************************
%* - Function of STAPMAT in initialization phase                   *
%*                                                                 *
%* - Purpose:                                                      *
%*     Read input file of STAPMAT                                  *
%*                                                                 *
%* - Call procedures:                                              *
%*     SRC/Initiation/ReadFile.m - InitBasicData()                 *
%*                                                                 *
%* - Called by :                                                   *
%*     stapmat.m                                                   *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.21                *
%*                                                                 *
%* *****************************************************************

function ReadFile(fname)
 % fname = strcat('.\Data\', fname);           % Deal the filename
fname = strcat('./Data/', fname);           % Deal the filename

% Get global class
global cdata;
global sdata;

% Open files
cdata.IIN = fopen(fname, 'r');

% Begin Read input file
fprintf('Input phase ...\n\n');

% the first time stamp
cdata.TIM = zeros(5, 6, 'double');
cdata.TIM(1,:) = clock;

IIN = cdata.IIN;
%% Read Control data
cdata.HED = fgetl(IIN);

tmp = str2num(fgetl(IIN));
cdata.NUMNP = int64(tmp(1));
cdata.NUMEG = int64(tmp(2));
cdata.NLCASE = int64(tmp(3));
cdata.MODEX = int64(tmp(4));%%如果为3，则采用显式的速度verlet格式求解。
cdata.numofeig = int64(tmp(5));%%增加了有关特征值求解的三个控制变量
cdata.blocksize = int64(tmp(6));
cdata.randominitial = int64(tmp(7));

if (cdata.NUMNP == 0) return; end
cdata.NPAR = zeros(3, cdata.NUMEG, 'int64');

%% Read Dynamics Parameters(Damping Matrix & Generalized Alpha Method & TimeStep等等等等)
sdata.MassType = 1;% default: use lumped mass
if(cdata.MODEX == 2||3)
    tmp = str2num(fgetl(IIN));
    sdata.ENDTIME = double(tmp(1));
    sdata.CparaM =  double(tmp(2));
    sdata.CparaK =  double(tmp(3));
    sdata.DyAlphaF = double(tmp(4));
    sdata.DyAlphaM = double(tmp(5));
    sdata.DyGamma = double(tmp(6));
    sdata.DyBeta = double(tmp(7));
    sdata.TimeStep = double(tmp(8));
    if(length(tmp)==9)
        if(tmp(9)==1||2)
            sdata.MassType = tmp(9);
        else
            error(' *** ERROR *** WRONG MassType');
        end
    end

    sdata.NSTEPS = ceil(sdata.ENDTIME/sdata.TimeStep)+1;
end

%% Read nodal point data
InitBasicData();
% Define local variables to speed
ID = sdata.ID; X = sdata.X; Y = sdata.Y; Z = sdata.Z;
for i = 1:cdata.NUMNP
    tmp = str2num(fgetl(IIN));
    ID(1, i) = int64(tmp(2));
    ID(2, i) = int64(tmp(3));
    ID(3, i) = int64(tmp(4));
    X(i) = double(tmp(5));
    Y(i) = double(tmp(6));
    Z(i) = double(tmp(7));
end
sdata.ID = ID; sdata.X = X; sdata.Y = Y; sdata.Z = Z;

%% Compute the number of equations
sdata.IDOrigin = ID;
NEQ = 0;
for N=1:cdata.NUMNP
    for I=1:3
        if (ID(I,N) == 0)
            NEQ = NEQ + 1;
            ID(I,N) = NEQ;
        else
            ID(I,N) = 0;
        end
    end
end
sdata.ID = ID;
sdata.NEQ = NEQ;

sdata.MASS = zeros(sdata.NEQ, 1, 'double');%%要质量阵的空间，要刚度阵的空间的代码在各单元的“trussstiff”中
sdata.lamda = zeros(1,cdata.numofeig,'double');%%这三行是要储存特征对的空间
sdata.omega = zeros(1,cdata.numofeig,'double');
sdata.PHI = zeros(sdata.NEQ, cdata.numofeig,'double');


%% Read Dynamic initial value (MODEX == 2或3)
if(cdata.MODEX == 2||3)
    tmp = str2num(fgetl(IIN));
    InitialDis = int64(tmp(1));
    InitialVol = int64(tmp(2));
    InitialAcc = int64(tmp(3));
    sdata.InitialDis = InitialDis;

    sdata.Ddis = zeros(sdata.NEQ, 1, 'double');
    sdata.Dvol = zeros(sdata.NEQ, 1, 'double');
    sdata.Dacc = zeros(sdata.NEQ, 1, 'double');

    if(InitialDis == 1)
        Ddis = sdata.Ddis;
        for node = 1:cdata.NUMNP
            tmp = str2num(fgetl(IIN));
            for dire = 1:3
                if(ID(dire,node))
                    Ddis(ID(dire,node)) = tmp(dire);
                end
            end
        end
        sdata.Ddis = Ddis;
    end

    if(InitialVol == 1)
        Dvol = sdata.Dvol;
        for node = 1:cdata.NUMNP
            tmp = str2num(fgetl(IIN));
            for dire = 1:3
                if(ID(dire,node))
                    Dvol(ID(dire,node)) = tmp(dire);
                end
            end
        end
        sdata.Dcol = Dvol;
    end

    if(InitialAcc == 1)
        Dacc = sdata.Dacc;
        for node = 1:cdata.NUMNP
            tmp = str2num(fgetl(IIN));
            for dire = 1:3
                if(ID(dire,node))
                    Dacc(ID(dire,node)) = tmp(dire);
                end
            end
        end
        sdata.Dacc = Dacc;
    end

end
%% Read load data
% Init control data
NLCASE = cdata.NLCASE;
sdata.R = zeros(NEQ, NLCASE, 'double');
R = sdata.R;

% Read data
for N = 1:cdata.NLCASE
    tmp = str2num(fgetl(IIN));
    if(isempty(cdata.LL)~=1 && cdata.LL>=int64(tmp(1)))
        error(' *** ERROR *** LOAD CASES ARE NOT IN ORDER');
    end
    cdata.LL = int64(tmp(1));
    cdata.NLOAD = int64(tmp(2));
    NLOAD = cdata.NLOAD;
%   Init load data
    sdata.NOD = zeros(NLOAD, 1, 'int64');
    sdata.IDIRN = zeros(NLOAD, 1, 'int64');
    sdata.FLOAD = zeros(NLOAD, 1, 'double');
    NOD = sdata.NOD; IDIRN = sdata.IDIRN; FLOAD = sdata.FLOAD;
    
%   Read load data
    for I = 1:NLOAD
        tmp = str2num(fgetl(IIN));
        NOD(I) = int64(tmp(1));
        IDIRN(I) = int64(tmp(2));
        FLOAD(I) = double(tmp(3));
    end
    if (cdata.MODEX == 0) return; end
    
%   Compute load vector
    for L = 1:NLOAD
        II = ID(IDIRN(L), NOD(L));
        if (II > 0) R(II, N) = R(II, N) + FLOAD(L); end
    end
    sdata.NOD = NOD; sdata.IDIRN = IDIRN; sdata.FLOAD = FLOAD; sdata.R = R;
end

%% Read Time-Variant Load data
if(cdata.MODEX == 2||3)
     tmp = str2num(fgetl(IIN));
     sdata.NVNL = int64(tmp(1));

     sdata.DLDC = zeros(1,sdata.NVNL,'int64');
     sdata.DNOD = zeros(1,sdata.NVNL,'int64');
     sdata.DDIRE = zeros(1,sdata.NVNL,'int64');
     DLDC = sdata.DLDC; DNOD = sdata.DNOD; DDIRE = sdata.DDIRE;
     sdata.DVNL = zeros(sdata.NSTEPS,sdata.NVNL,'double');
     DVNL = sdata.DVNL;

     for I = 1:sdata.NVNL
         tmp = str2num(fgetl(IIN));
         DLDC(I) = int64(tmp(1));
         DNOD(I) = int64(tmp(2));
         DDIRE(I) = int64(tmp(3));

         tmp = str2num(fgetl(IIN));
         if (length(tmp) == sdata.NSTEPS)
             for J = 1:sdata.NSTEPS
                 DVNL(J,I) = tmp(J);
             end
         else
             error('Time-Varient Load Input DOESN''T MATCH the Time Step!');
         end
     end

     sdata.DLDC = DLDC; sdata.DNOD = DNOD; sdata.DDIRE = DDIRE; sdata.DVNL = DVNL;

end

end



%% Functions


% InitBasicData
function InitBasicData()
global cdata;
global sdata;

% cdata.NPAR = zeros(10, 1, 'int64');

sdata.ID = zeros(3,cdata.NUMNP, 'int64');
sdata.X = zeros(cdata.NUMNP, 1, 'double');
sdata.Y = zeros(cdata.NUMNP, 1, 'double');
sdata.Z = zeros(cdata.NUMNP, 1, 'double');
end