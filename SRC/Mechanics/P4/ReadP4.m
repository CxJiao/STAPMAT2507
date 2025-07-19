%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     Read the element information of P4 element                  *
%*                                                                 *
%* - Call procedures:                                              *
%*     ReadP4.m - ReadMaterial()                                   *
%*     ReadP4.m - ReadElements()                                   *
%*                                                                 *
%* - Called by :                                                   *
%*     ./P4Stiff.m                                                 *
%*                                                                 *
%*                                                                 *
%* *****************************************************************

function ReadP4(NG)

% Read Material information
ReadMaterial(NG)

% Read Element information
ReadElements(NG)

% Read Pressure load for Plate element
ReadPressure()

% the second time stamp
global cdata;
cdata.TIM(2,:) = clock;

end

function ReadMaterial(NG)
global cdata;
global sdata;
% Get file pointers
IIN = cdata.IIN;
IOUT = cdata.IOUT;

if (cdata.NPAR(3,NG) == 0) 
    cdata.NPAR(3,NG) = 1; 
end
fprintf(IOUT, '\n M A T E R I A L   D E F I N I T I O N\n');
fprintf(IOUT, '\n NUMBER OF DIFFERENT SETS OF MATERIAL\n');
fprintf(IOUT, ' AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . = %10d\n', cdata.NPAR(3,NG));
fprintf(IOUT, '     SET           YOUNG''S       POISSON        PLATE          DENSITY\n');
fprintf(IOUT, ' NUMBER     MODULUS        RATIO       THICKNESS\n');
fprintf(IOUT, '                            E                   nu                   t                   rho\n');


% Read material datas
sdata.NUME = cdata.NPAR(2,NG);
sdata.NUMMAT = cdata.NPAR(3,NG);
NUMMAT = cdata.NPAR(3,NG);
sdata.E = zeros(NUMMAT, 1, 'double');
sdata.t = zeros(NUMMAT, 1, 'double');
sdata.nu = zeros(NUMMAT, 1, 'double');
sdata.rho = zeros(NUMMAT, 1, 'double');
for I = 1:cdata.NPAR(3,NG)
    tmp = str2num(fgetl(IIN));
    N = round(tmp(1));
    sdata.E(N) = tmp(2);
    sdata.nu(N) = tmp(3);
    sdata.t(N) = tmp(4); 
    sdata.rho(N) = tmp(5);
    fprintf(IOUT, '    %5d          %9.2e       %9.2e      %9.2e       %9.2e\n', ...
        N, tmp(2), tmp(3), tmp(4), tmp(5));
end

end

% Read elements information
function ReadElements(NG)
global cdata;
global sdata;

% Get file pointer
IIN = cdata.IIN;
IOUT = cdata.IOUT;

fprintf(IOUT, '\n\n E L E M E N T   I N F O R M A T I O N\n');
fprintf(IOUT, '\n    ELEMENT                 NODE             MATERIAL      PRESSURE\n');
fprintf(IOUT, '  NUMBER-N        A     B     C     D          SET                  SET\n');

% Get Position data
NUME = cdata.NPAR(2,NG);
sdata.XYZ = zeros(12, NUME, 'double');
sdata.ELNOD = zeros(sdata.NNODE,NUME,'int64');
sdata.MATP = zeros(1, NUME, 'int64');
sdata.LM = zeros(12, NUME, 'double');
sdata.MHT = zeros(sdata.NEQ, 1, 'int64');
sdata.PRSID = zeros(1, NUME, 'double');
X = sdata.X; Y = sdata.Y; Z = sdata.Z; ID = sdata.ID;
MATP = sdata.MATP; PRSID=sdata.PRSID;
P4NODE=zeros(4,1);
XYZ=zeros(12,1,'double');
LM=zeros(12,1,'double');
for N = 1:NUME
    tmp = str2num(fgetl(IIN));
    for i = 2:5  
        P4NODE(i-1) = round(tmp(i));  
    end
    MATP(N) = round(tmp(6)); 
    PRSID(N) = round(tmp(7));
    sdata.ELNOD(:,N)=P4NODE;
    %   Save element information
    for i = 1:4
        XYZ(3*i-2) = X(P4NODE(i));
        XYZ(3*i-1) = Y(P4NODE(i));
        XYZ(3*i) = Z(P4NODE(i));
    end
    sdata.XYZ(:,N) = XYZ;
    fprintf(IOUT,'%10d              ',N);
    fprintf(IOUT,'%6d',P4NODE);
    fprintf(IOUT,'           %2d            ',MATP(N));
    fprintf(IOUT,'           %2d            \n',PRSID(N));

%   Compute connectivity matrix
    for i = 1:4
        LM(3*i-2)=ID(1, P4NODE(i));
        LM(3*i-1)=ID(2, P4NODE(i));
        LM(3*i)=ID(3, P4NODE(i));
    end

%   Updata column heights and bandwidth
    ColHt(LM)
    sdata.LM(:,N) = LM;
end
sdata.MATP = MATP;
sdata.PRSID = PRSID;

end

function ReadPressure()
global cdata;
global sdata;

IIN = cdata.IIN;
IOUT = cdata.IOUT;

PRSN = max(sdata.PRSID);
NSTEPS=sdata.NSTEPS;
sdata.PRESS = zeros(NSTEPS, PRSN, 'double');
PRESS=sdata.PRESS;

if PRSN > 0
    for i = 1:PRSN
        tmp = str2num(fgetl(IIN));
        if (length(tmp) == NSTEPS)
            for j = 1:NSTEPS
                PRESS(j,i) = tmp(j);
            end
        else
            error('Time-Varient Pressure DOESN''T MATCH the Time Step!');
        end
    end
end

sdata.PRESS = PRESS;
fprintf(IOUT, '\n\n Time-Variant pressure load %d types have been read\n',i);

end
