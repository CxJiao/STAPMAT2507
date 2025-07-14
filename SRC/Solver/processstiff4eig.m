%* *****************************************************************
%* - Function of STAPMAT in Solver phase                           *
%*                                                                 *
%* - Purpose:                                                      *
%*     输入组装之后的总体刚度阵stiff4eig与质量阵MASS        *
%*     进行MASS^(-1/2）*stiff4eig*MASS^(-1/2）后再LDLT分解，使问题变为一个标准特征值问题，更改后的量仍存储在stiff4eig中   *
%* - Call procedures:（暂时没有更新）                                              *
%*     ./LDLTFactor.m            - LDLTFactor()                    *
%*     Solve.m                   - Stiff2Sparse()                  *
%*     ./ColSol.m                - ColSol()                        *  
%*     Solve.m                   - WriteDis()                      *
%*     SRC/Mechanics/GetStress.m - GetStress()                     *
%*                                                                 *
%* - Called by :     （暂时没有更新）                                              *
%*     stapmat.m                                                   *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.22                *
%*                                                                 *
%* *****************************************************************

function processstiff4eig()

global sdata;

A = sdata.stiff4eig; 
M = sdata.MASS;
MAXA = sdata.MAXA; NEQ = sdata.NEQ; 

if min(M)<=0
    error(['STOP 质量阵的对角线出现零元或负元']);
end
M=M.^(-1/2);


%%一列一列地变换A的元素
for N = 1:NEQ %N代表矩阵的列，即某元素位于第N列
    KN = MAXA(N);
    KU = MAXA(N + 1);
    total=KU-KN;%本列总共非0，即需要更改的元素总数
    for NN=1:total
        %该元素在A中的地址为MAXA(N)+NN-1
        %该元素位于的行是N+1-NN
        A(MAXA(N)+NN-1)=A(MAXA(N)+NN-1)*M(N)*M(N+1-NN);%%MASS^(-1/2）*stiff4eig*MASS^(-1/2）
    end
end

%照搬LDLT分解的代码，对A进行LDLT分解，方便逆迭代求标准大型特征值问题
for N = 1:NEQ
    KN = MAXA(N);
    KL = KN + 1;
    KU = MAXA(N + 1) - 1;
    KH = KU - KL;
    
    if (KH > 0)
        K = N - KH;
        IC = 0;
        KLT = KU;
        for J = 1:KH
            IC = IC + 1;
            KLT = KLT - 1;
            KI = MAXA(K);
            ND = MAXA(K+1) - KI - 1;
            if (ND > 0)
                KK = min(IC, ND);
                C = 0.0;
                for L = 1:KK C = C+A(KI+L)*A(KLT+L); end
                A(KLT) = A(KLT) - C;
            end
            K = K + 1;
        end
    end
    
    if (KH >= 0)
        K = N;
        B = 0.0;
        for KK = KL:KU
            K = K - 1;
            KI = MAXA(K);
            C = A(KK) / A(KI);
            B = B + C*A(KK);
            A(KK) = C;
        end
        A(KN) = A(KN) - B;
    end
    
    if (A(KN) <= 0)
        error(['STOP - Stiffness matrix is not positive definite\n' ...
            'Nonpositive number for equation %8d is %20.12e\n'], N, A(KN));
    end
end

sdata.stiff4eig = A;%%将处理完的A返回到全局变量中，后续只需要处理stiff4eig的标准特征值问题

end


