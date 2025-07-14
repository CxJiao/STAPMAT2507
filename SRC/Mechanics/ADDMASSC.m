%* *****************************************************************
%* - Function of STAPMAT in compatible phase                       *
%*                                                                 *
%* - Purpose:                                                      *
%*     To assemble mass into global concentrated mass              *
%*                                                                 *
%* - Call procedures: None                                         *
%*                                                                 *
%* - Called by :                                                   *
%*     ./Truss/TrussStiff.m - Assemble()                           *
%*     ./H20/H20Stiff.m - Assemble()                               *
%*                                                                 *
%* *****************************************************************

function ADDMASSC(M, LM)

% Get global data
global sdata;
MASSC = sdata.MASSC;
MAXA = sdata.MAXA;
ND = sdata.NDOF * sdata.NNODE;
for J = 1:ND
    JJ = LM(J);
    if (JJ > 0)
        for I = 1:J
            II = LM(I);
            if (II > 0)
                if (JJ > II) KK = MAXA(JJ) + JJ - II;
                else KK = MAXA(II) + II - JJ; end
                MASSC(KK) = MASSC(KK) + M(I, J);
            end
        end
    end
end
sdata.MASSC = MASSC;

end