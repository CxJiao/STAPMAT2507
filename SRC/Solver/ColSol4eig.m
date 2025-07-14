%* *****************************************************************
%* - Function of STAPMAT in Solver phase                           *
%*                                                                 *
%* - Purpose:                                                      *
%*     已知一个已经QR分解好的矩阵stiff4eig，输入矩阵X，求A\X       *
%*     专门用于逆特征值求解的解方程部分                                                            *
%* - Call procedures: None       (暂时没有更改）                                  *
%*                                                                 *
%* - Called by :                                                   *
%*     ./Solve.m                                                   *
%*                                                                 *
%* - Programmed in Fortran 90 by Xiong Zhang                       *
%*                                                                 *
%* - Adapted to Matlab by:                                         *
%*     Yan Liu, Computational Dynamics Group, School of Aerospace  *
%*     Engineering, Tsinghua University, 2019.02.22                *
%*                                                                 *
%* *****************************************************************

function [Xout]=ColSol4eig(X)

% Get global data
global sdata;
A = sdata.stiff4eig; MAXA = sdata.MAXA; 
NEQ = sdata.NEQ; NWK = sdata.NWK; 
NNM = NEQ + 1;

Xout=zeros(NEQ,size(X,2));
for NUM=1:size(X,2)

    RR = X(:,NUM);

    % Reduce right-hand-side load vector
    for N = 1:NEQ
        KL = MAXA(N) + 1;
        KU = MAXA(N+1) - 1;
        if (KU-KL >= 0)
            K = N;
            C = 0.0;
            for KK = KL:KU
                K = K - 1;
                C = C + A(KK) * RR(K);
            end
            RR(N) = RR(N) - C;
        end
    end

    % Back-Substitute
    for N = 1:NEQ
        K = MAXA(N);
        RR(N) = RR(N) / A(K);
    end

    if (NEQ == 1) 
        return; 
    end

    N = NEQ;
    for L = 2:NEQ
        KL = MAXA(N) + 1;
        KU = MAXA(N+1) - 1;
        if (KU-KL >= 0)
            K = N;
            for KK = KL:KU
                K = K - 1;
                RR(K) = RR(K) - A(KK)*RR(N);
            end
        end
        N = N - 1;
    end
Xout(:,NUM)=RR;

end
end