%* *****************************************************************
%* - Function of STAPMAT in lumped mass phase                      *
%*                                                                 *
%* - Purpose:                                                      *
%*     To assemble pressure node load into global load vector DVNL *
%*                                                                 *
%* - Call procedures: None                                         *
%*                                                                 *
%* - Called by :                                                   *
%*     ./P4/P4Stiff.m - Assemble()                                 *
%*                                                                 *
%* *****************************************************************

function ADDPRESS(PE, N)
% Get global data
global sdata;
global cdata;
PRESS = sdata.PRESS;
PRSID = sdata.PRSID;
ID = sdata.ID;
LMN = sdata.LM(:,N);
if cdata.MODEX == 1
    % only the pressure at time = 0 is valid for static analysis
    R = sdata.R;
    PTYPE=PRSID(N);
    if PTYPE > 0
        LOAD=PRESS(1,PTYPE);
        for j=1:12
            if LMN(j) > 0
                R(LMN(j)) = R(LMN(j)) + LOAD*PE(j);
            end
        end
    end
    sdata.R = R;
elseif cdata.MODEX == 2||3
    if sdata.NVNL==0
        % initial the dynamic parameter
        sdata.DLDC = zeros(1,0,'int64');
        sdata.DNOD = zeros(1,0,'int64');
        sdata.DDIRE = zeros(1,0,'int64');
        sdata.DVNL = zeros(sdata.NSTEPS,0,'double');
    end
    NVNL = sdata.NVNL;
    DLDC = sdata.DLDC;
    DNOD = sdata.DNOD;
    DDIRE = sdata.DDIRE;
    DVNL = sdata.DVNL;
    
    PTYPE=PRSID(N);
    if PTYPE > 0
        LOADS=PRESS(:,PTYPE);
        for j=1:12
            if LMN(j) > 0
                [idx,idy] = find(ID == LMN(j));
                DNOD = [DNOD idy];
                DDIRE = [DDIRE idx];
                DVNL = [DVNL PE(j)*LOADS];
                NVNL = NVNL + 1;
            end
        end
    end
    sdata.DLDC = ones(1,NVNL); % 这里偷个懒，全部都填工况1
    sdata.DVNL = DVNL;
    sdata.NVNL = NVNL;
    sdata.DNOD = DNOD;
    sdata.DDIRE = DDIRE;
end

end
