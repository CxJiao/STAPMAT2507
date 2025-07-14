%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     To assemble element stiffness into global stiffness         *
%*                                                                 *
%* - Call procedures: None                                         *
%*                                                                 *
%* - Called by :                                                   *
%*     ./Truss/TrussStiff.m - Assemble()                           *
%*     ./H20/H20Stiff.m - Assemble()                               *
%*                                                                 *
%* - Programmed in Fortran 90 by Xiong Zhang                       *
%*                                                                 *
%* - Adapted to Matlab by:                                         *
%*     LeiYang Zhao, Yan Liu, Computational Dynamics Group,        *
%*     School of Aerospace Engineering, Tsinghua University,       *
%*     2019.02.22                                                  *
%*                                                                 *
%* *****************************************************************

function ADDBAN(S, LM)

% Get global data
global sdata;
MAXA = sdata.MAXA; 
STIFF = sdata.STIFF;
ND = sdata.NDOF * sdata.NNODE;
for J = 1:ND
    JJ = LM(J);
    if (JJ > 0)
        for I = 1:J
            II = LM(I);
            if (II > 0)
                if (JJ > II) KK = MAXA(JJ) + JJ - II;
                else KK = MAXA(II) + II - JJ; end
                STIFF(KK) = STIFF(KK) + S(I, J);
            end
        end
    end
end

sdata.STIFF = STIFF;
sdata.stiff4eig = STIFF;%这个刚度阵将来专门做特征值求解用

end