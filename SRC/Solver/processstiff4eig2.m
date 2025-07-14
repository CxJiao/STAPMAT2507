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

function processstiff4eig2()

global sdata;

A = sdata.STIFF; %%已经LDLT分解好的矩阵
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
        if NN==1
            A(MAXA(N)+NN-1)=A(MAXA(N)+NN-1)*M(N)*M(N);%%LDLT的对角元直接乘M就行
        else
            A(MAXA(N)+NN-1)=A(MAXA(N)+NN-1)*M(N)/M(N+1-NN);%%非对角元这样处理
        end
    end
end


sdata.stiff4eig = A;%%将处理完的A返回到全局变量中，后续只需要处理stiff4eig的标准特征值问题

end


