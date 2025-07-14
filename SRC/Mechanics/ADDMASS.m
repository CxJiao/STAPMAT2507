%* *****************************************************************
%* - Function of STAPMAT in lumped mass phase                      *
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

function ADDMASS(M, LM)

% Get global data
global sdata;
MASS = sdata.MASS;
ND = sdata.NDOF * sdata.NNODE;
for J = 1:ND
    JJ = LM(J);
    if (JJ > 0)
        MASS(JJ) = MASS(JJ) + M(J);
    end
end

sdata.MASS = MASS;

end