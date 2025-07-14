%* *****************************************************************
%* - Function of STAPMAT in solver phase                           *
%*                                                                 *
%* - Purpose:                                                      *
%*     To calculate stresses                                       *
%*                                                                 *
%* - Call procedures: None                                         *
%*                                                                 *
%* - Called by :                                                   *
%*     SRC/Solver/GetStress.m                                      *
%*                                                                 *
%* *****************************************************************

function H20Stress(NUM,NG)

% Get global data
global cdata;
global sdata;

IOUT = cdata.IOUT;
NUME = sdata.NUME2; MATP = sdata.MATP2; 
E = sdata.E2; nu = sdata.nu; LM = sdata.LM2;
U = sdata.DIS(:, NUM); XYZ = sdata.XYZ2; GC3 = sdata.GC3;

xyz_int = zeros(27*NUME,3);
stress = zeros(27*NUME,6);
index = 1;

fprintf(IOUT, ['\n\n  S T R A I N  C A L C U L A T I O N S  F O R  ' ...
    'E L E M E N T  G R O U P %4d\n\n' ... 
    '                LOCATION                             x               y' ...
    '               z              xy              yz              zx\n'], NG);

for N = 1:NUME
    % MTYPE = MATP(N); v = nu(N);
    MTYPE = MATP(N); v = nu(MTYPE);

    D = E(MTYPE) / (1 + v) / (1 - 2*v) * [1-v v v 0 0 0;...
        v 1-v v 0 0 0;v v 1-v 0 0 0;0 0 0 (1-2*v)/2 0 0;...
        0 0 0 0 (1-2*v)/2 0;0 0 0 0 0 (1-2*v)/2];

    % calculate displacements of element nodes
    dis = zeros(60,1);
    for i = 1:60
        if(LM(i,N)>0)
            dis(i) = U(LM(i,N));
        end
    end

    xyz = XYZ(:, N);
    
    for i = 1:3
        for j = 1:3
            for k = 1:3
                xyz_int(index, :) = N_cal(GC3(i),GC3(j),GC3(k)) * xyz;% integration point
               
                strain = B_cal(GC3(i),GC3(j),GC3(k),N) * dis;% strain
                stress(index, :) = D * strain;

                fprintf(IOUT, '%5e  ',xyz_int(index, :));
                fprintf(IOUT, '   %13.6e',strain);
                fprintf(IOUT, '\n');

                index = index + 1;
            end
        end
    end

end


fprintf(IOUT, ['\n\n  S T R E S S  C A L C U L A T I O N S  F O R  ' ...
    'E L E M E N T  G R O U P %4d\n\n' ... 
    '                LOCATION                             x               y' ...
    '               z              xy              yz              zx\n'], NG);

for index = 1:size(stress,1)
    fprintf(IOUT, '%5e  ',xyz_int(index, :));
    fprintf(IOUT, '   %13.6e',stress(index, :));
    fprintf(IOUT, '\n');
end


end