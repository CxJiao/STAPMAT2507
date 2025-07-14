%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     Read the element information of H20                         *
%*                                                                 *
%* - Call procedures:                                              *
%*     ReadH20.m - ReadMaterial()                                  *
%*     ReadH20.m - ReadElements() - ReadAssembleLoads()            *
%*                                                                 *
%* - Called by :                                                   *
%*     ./H20Stiff.m                                                *
%*                                                                 *
%* *****************************************************************

function ReadH20(NG)

% Read Material information
ReadMaterial(NG)

% Read Element information
ReadElements(NG)

% the second time stamp
global cdata;
cdata.TIM(2,:) = clock;

end


% ----------------------- Functions -----------------------------------
% Read Material information
function ReadMaterial(NG)

global cdata;
global sdata;
% Get file pointers
IIN = cdata.IIN;
IOUT = cdata.IOUT;

if (cdata.NPAR(3,NG) == 0) cdata.NPAR(3,NG) = 1; end
fprintf(IOUT, '\n M A T E R I A L   D E F I N I T I O N\n');
fprintf(IOUT, '\n NUMBER OF DIFFERENT SETS OF MATERIAL. .( NPAR(3) ) . . = %10d\n' ,...
    cdata.NPAR(3,NG));
fprintf(IOUT, '  SET       YOUNG''S        POISSON     DENSITY\n');
fprintf(IOUT, ' NUMBER     MODULUS         RATIO\n');
fprintf(IOUT, '               E             nu           rho\n');


% Read material datas
sdata.NUME2 = cdata.NPAR(2,NG);
sdata.NUMMAT = cdata.NPAR(3,NG);
NUMMAT = cdata.NPAR(3,NG);
sdata.E2 = zeros(NUMMAT, 1, 'double');
sdata.nu = zeros(NUMMAT, 1, 'double');
sdata.rho = zeros(NUMMAT, 1, 'double');
for I = 1:cdata.NPAR(3,NG)
    tmp = str2num(fgetl(IIN));
    N = round(tmp(1));
    sdata.E2(N) = tmp(2);
    sdata.nu(N) = tmp(3);
    sdata.rho(N) = tmp(4);
    fprintf(IOUT, '%5d    %12.5e  %14.6e  %14.6e\n', N, tmp(2), tmp(3), tmp(4));
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
fprintf(IOUT, ['\n      ELEMENT     NODE                                 ' ...
    '                                              MATERIAL\n']);
fprintf(IOUT, ['      NUMBER-N       A     B     C     D     E     F     G     H     I     J' ...
    '     K     L     M     N     O     P     Q     R     S     T    SET NUMBER\n']);


% Get Position data
NUME = cdata.NPAR(2,NG);
sdata.XYZ2 = zeros(60, NUME, 'double');
sdata.MATP2 = zeros(NUME, 1, 'int64');                 % the type of material
sdata.LM2 = zeros(60, NUME, 'double');                  % connectivity matrix
sdata.MHT = zeros(sdata.NEQ, 1, 'int64');
X = sdata.X; Y = sdata.Y; Z = sdata.Z; ID = sdata.ID;
XYZ = sdata.XYZ2; MATP = sdata.MATP2; LM = sdata.LM2;
H20NODE=zeros(20,1);

for N = 1:NUME
    tmp = str2num(fgetl(IIN));
    for i = 2:21  H20NODE(i-1) = round(tmp(i));  end
    MTYPE = round(tmp(22)); 
    NF = round(tmp(23)); NV = round(tmp(24));

    %   Save element information
    for i = 1:20
        XYZ(3*i-2,N) = X(H20NODE(i));
        XYZ(3*i-1,N) = Y(H20NODE(i));
        XYZ(3*i,N) = Z(H20NODE(i));
    end
    MATP(N) = MTYPE;
    sdata.XYZ2 = XYZ;

    fprintf(IOUT,'%10d      ',N);
    fprintf(IOUT,'%6d',H20NODE);
    fprintf(IOUT,'     %5d\n',MTYPE);

%   Compute connectivity matrix
    for i = 1:20
        LM(3*i-2, N)=ID(1, H20NODE(i));
        LM(3*i-1, N)=ID(2, H20NODE(i));
        LM(3*i, N)=ID(3, H20NODE(i));
    end

%   Updata column heights and bandwidth
    ColHt(LM(:, N))

    sdata.LM2 = LM;
    ReadAssembleLoads(NF,NV,N,H20NODE);
end
sdata.MATP2 = MATP;

% Clear the memory of X, Y, Z
if( NG == cdata.NUMEG )
    sdata.X = double(0);
    sdata.Y = double(0);
    sdata.Z = double(0);
end

end



% Read Load information
function ReadAssembleLoads(NF,NV,NG,NODE)
% NF, NV: amounts of face loads and volume loads on this element
% NG: element index
% XYZ: coordinates of this element nodes

global cdata;
global sdata;

% Get file pointer
IIN = cdata.IIN;
GC3 = sdata.GC3; GW3 = sdata.GW3;
R = sdata.R; LM = sdata.LM2;

for i = 1:NF+NV
    tmp = str2num(fgetl(IIN));
    LoadCase = int64(tmp(1));
    LoadType = int64(tmp(2));
    Rx = zeros(60,1);
    
    if(LoadType==1)% Face Load
        FaceOn = [int64(tmp(3)) int64(tmp(4)) int64(tmp(5)) int64(tmp(6))];% Global Coordinate
        faceOn = zeros(1,4);% Local Coordinate
        FIDIRN = int64(tmp(7));
        FMAG = double(tmp(8));

        FFLOAD = zeros(3,1);
        if( FIDIRN==1 || FIDIRN==2 || FIDIRN==3 )
            FFLOAD(FIDIRN) = FMAG;
        else
            error(' *** ERROR *** No Such direction Of Load');
        end

        for loc = 1:4
            faceOn(loc) = find(NODE == FaceOn(loc));
        end
        faceOn = sort(faceOn);

        if(all(faceOn == [1 4 5 8]))
            for j = 1:3
                for k = 1:3
                    JJ = JACOBI(-1,GC3(j),GC3(k),NG);
                    dS = sqrt((JJ(2,2)*JJ(3,3)-JJ(2,3)*JJ(3,2))^2+(JJ(2,3)*JJ(3,1)-JJ(3,3)*JJ(2,1))^2+...
                        (JJ(2,1)*JJ(3,2)-JJ(2,2)*JJ(3,1))^2);
                    Rx = Rx + GW3(j)*GW3(k)*(N_cal(-1,GC3(j),GC3(k)))'*FFLOAD*dS;
                    % Rx = Rx + GW3(j)*GW3(k)*(N_cal(-1,GC3(j),GC3(k)))'*det(JJ(2:3,2:3))*FFLOAD;
                end
            end
        elseif(all(faceOn == [2 3 6 7]))
            for j = 1:3
                for k = 1:3
                    JJ = JACOBI(1,GC3(j),GC3(k),NG);
                    dS = sqrt((JJ(2,2)*JJ(3,3)-JJ(2,3)*JJ(3,2))^2+(JJ(2,3)*JJ(3,1)-JJ(3,3)*JJ(2,1))^2+...
                        (JJ(2,1)*JJ(3,2)-JJ(2,2)*JJ(3,1))^2);
                    Rx = Rx + GW3(j)*GW3(k)*(N_cal(1,GC3(j),GC3(k)))'*FFLOAD*dS;
                    % Rx = Rx + GW3(j)*GW3(k)*(N_cal(1,GC3(j),GC3(k)))'*det(JJ(2:3,2:3))*FFLOAD;
                end
            end
        elseif(all(faceOn == [1 2 5 6]))
            for j = 1:3
                for k = 1:3
                    JJ = JACOBI(GC3(j),-1,GC3(k),NG);
                    dS = sqrt((JJ(1,1)*JJ(3,3)-JJ(1,3)*JJ(3,1))^2+(JJ(1,1)*JJ(3,2)-JJ(1,2)*JJ(3,1))^2+...
                        (JJ(1,2)*JJ(3,3)-JJ(1,3)*JJ(3,2))^2);
                    Rx = Rx + GW3(j)*GW3(k)*(N_cal(GC3(j),-1,GC3(k)))'*FFLOAD*dS;
                    % Rx = Rx + GW3(j)*GW3(k)*(N_cal(GC3(j),-1,GC3(k)))'*det(JJ([1,3],[1,3]))*FFLOAD;
                end
            end
        elseif(all(faceOn == [3 4 7 8]))
            for j = 1:3
                for k = 1:3
                    JJ = JACOBI(GC3(j),1,GC3(k),NG);
                    dS = sqrt((JJ(1,1)*JJ(3,3)-JJ(1,3)*JJ(3,1))^2+(JJ(1,1)*JJ(3,2)-JJ(1,2)*JJ(3,1))^2+...
                        (JJ(1,2)*JJ(3,3)-JJ(1,3)*JJ(3,2))^2);
                    Rx = Rx + GW3(j)*GW3(k)*(N_cal(GC3(j),1,GC3(k)))'*FFLOAD*dS;
                    % Rx = Rx + GW3(j)*GW3(k)*(N_cal(GC3(j),1,GC3(k)))'*det(JJ([1,3],[1,3]))*FFLOAD;
                end
            end
        elseif(all(faceOn == [1 2 3 4]))
            for j = 1:3
                for k = 1:3
                    JJ = JACOBI(GC3(j),GC3(k),-1,NG);
                    dS = sqrt((JJ(1,1)*JJ(2,2)-JJ(1,2)*JJ(2,1))^2+(JJ(1,1)*JJ(2,3)-JJ(1,3)*JJ(2,1))^2+...
                        (JJ(1,2)*JJ(2,3)-JJ(2,2)*JJ(1,3))^2);
                    Rx = Rx + GW3(j)*GW3(k)*(N_cal(GC3(j),GC3(k),-1))'*FFLOAD*dS;
                    % Rx = Rx + GW3(j)*GW3(k)*(N_cal(GC3(j),GC3(k),-1))'*det(JJ(1:2,1:2))*FFLOAD;
                end
            end
        elseif(all(faceOn == [5 6 7 8]))
            for j = 1:3
                for k = 1:3
                    JJ = JACOBI(GC3(j),GC3(k),1,NG);
                    dS = sqrt((JJ(1,1)*JJ(2,2)-JJ(1,2)*JJ(2,1))^2+(JJ(1,1)*JJ(2,3)-JJ(1,3)*JJ(2,1))^2+...
                        (JJ(1,2)*JJ(2,3)-JJ(2,2)*JJ(1,3))^2);
                    Rx = Rx + GW3(j)*GW3(k)*(N_cal(GC3(j),GC3(k),1))'*FFLOAD*dS;
                    % Rx = Rx + GW3(j)*GW3(k)*(N_cal(GC3(j),GC3(k),1))'*det(JJ(1:2,1:2))*FFLOAD;
                end
            end
        else
            error(' *** ERROR *** Wrong Face Load Input');
        end

        
    elseif(LoadType==2)% Volume Load
        VIDIRN = int64(tmp(3));% coordinate axis(direction)
        VMAG = double(tmp(4));% magnitude

        VVLOAD = zeros(3,1);
        if( VIDIRN==1 || VIDIRN==2 || VIDIRN==3 )
            VVLOAD(VIDIRN) = VMAG;
        else
            error(' *** ERROR *** No Such direction Of Load');
        end

        for j = 1:3
            for k = 1:3
                for l = 1:3
                    Rx = Rx+GW3(j)*GW3(k)*GW3(l)*N_cal(GC3(j),GC3(k),GC3(l))'...
                        *det(JACOBI(GC3(j),GC3(k),GC3(l),NG))*VVLOAD;
                end
            end
        end

    else
        error(' *** ERROR *** No Such Type Of Load')
    end
    for p=1:60
        if(LM(p,NG) > 0) R(LM(p,NG),LoadCase) = R(LM(p,NG),LoadCase)+Rx(p); end
    end
end

sdata.R = R;

end