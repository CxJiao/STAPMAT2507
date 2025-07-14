%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     (for 1st element group) To calculate addresses of diagonal  *
%*     elements in banded matrix whose column heights are known    *
%*     (for else) And update global stiffness matrix               *                      
%*                                                                 *
%* - Call procedures: None                                         *
%*                                                                 *
%* - Called by :                                                   *
%*     ./Truss/TrussStiff.m                                        *
%*     ./H20/H20Stiff.m                                            *
%*                                                                 *
%* - Programmed in Fortran 90 by Xiong Zhang                       *
%*                                                                 *
%* - Adapted to Matlab by:                                         *
%*     LeiYang Zhao, Yan Liu, Computational Dynamics Group,        *
%*     School of Aerospace Engineering, Tsinghua University,       *
%*     2019.02.22                                                  *
%*                                                                 *
%* *****************************************************************

function Addres(NG)

% Get global data
global sdata;
global cdata;

NEQ = sdata.NEQ; MHT = sdata.MHT;
if( NG == 1 )
    sdata.MAXA = zeros(NEQ+1, 1, 'int64');
    MAXA = sdata.MAXA;

    MAXA(1) = 1;
    MAXA(2) = 2;
    MK = 0;

    if (NEQ > 1)
        for I = 2:NEQ
            if (MHT(I) > MK) MK = MHT(I); end
            MAXA(I+1) = MAXA(I) + MHT(I) + 1;
        end
    end

    sdata.MK = MK + 1;
    sdata.NWK = MAXA(NEQ+1) - MAXA(1);
    sdata.MAXA = MAXA;

else
    MAXA0 = sdata.MAXA;
    MHT0 = zeros(NEQ, 1, 'int64');
    MHT0(1) = 0;
    if (NEQ > 1)
        for I = 2:NEQ
            MHT0(I) = MAXA0(I+1)-MAXA0(I)-1;
        end
    end
    MHT1 = max(MHT0,MHT); sdata.MHT = MHT1;

    MK = 0; MAXA1 = zeros(NEQ+1, 1, 'int64');
    MAXA1(1) = 1; MAXA1(2) = 2;
    if (NEQ > 1)
        for I = 2:NEQ
            if (MHT1(I) > MK) MK = MHT1(I); end
            MAXA1(I+1) = MAXA1(I) + MHT1(I) + 1;
        end
    end
    sdata.MAXA = MAXA1; sdata.MK = MK+1;
    sdata.NWK = MAXA1(NEQ+1) - MAXA1(1);

    % update STIFF
    STIFF0 = sdata.STIFF;% former
    STIFF1 = zeros(sdata.NWK, 1, 'double');% updated
    STIFF1(1) = STIFF0(1);
    for I = 2:NEQ
        STIFF1(MAXA1(I):MAXA1(I)+MHT0(I)) = STIFF0(MAXA0(I):MAXA0(I+1)-1);
    end
    sdata.STIFF = STIFF1;

    if(sdata.MassType == 2)
        % update Compatible MASS
        MASS0 = sdata.MASSC;% former
        MASS1 = zeros(sdata.NWK, 1, 'double');% updated
        MASS1(1) = MASS0(1);
        for I = 2:NEQ
            MASS1(MAXA1(I):MAXA1(I)+MHT0(I)) = MASS0(MAXA0(I):MAXA0(I+1)-1);
        end
        sdata.MASSC = MASS1;
    end

end


% Write total system data
MM = round(sdata.NWK / NEQ);
fprintf(cdata.IOUT, ['\n\n  TOTAL SYSTEM DATA\n\n' ...
    '     NUMBER OF EQUATIONS . . . . . . . . . . . . . .(NEQ) = %10d\n' ...
    '     NUMBER OF MATRIX ELEMENTS . . . . . . . . . . .(NWK) = %10d\n' ...
    '     MAXIMUM HALF BANDWIDTH  . . . . . . . . . . . .(MK ) = %10d\n' ...
    '     MEAN HALF BANDWIDTH . . . . . . . . . . . . . .(MM ) = %10d\n'], ...
    NEQ, sdata.NWK, sdata.MK, MM);

end