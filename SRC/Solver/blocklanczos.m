%* *****************************************************************
%* - Function of STAPMAT in Solver phase                           *
%*                                                                 *
%* - Purpose:                                                      *
%*     输入已经转化为标准特征值问题且进行了LDLT分解的矩阵A        *
%*     要通过块lanczos变换法（块大小bsize）求出系统的前numeig个特征对，并将特征对保存进sdata中*
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

function blocklanczos()

global cdata;
global sdata;

reorjudge=1;%是否重正交化的开关，为1时在迭代过程中进行重正交化，为0是不进行重正交化，以后可以做一个主程序的接口
emsil=10^-6;%重正交化的系数

A = sdata.stiff4eig; %A已经是标准特征值问题的矩阵，且进行了LDLT分解，便于解方程逆迭代
M = sdata.MASS;
MAXA = sdata.MAXA; 
NEQ = sdata.NEQ; 

numeig=cdata.numofeig;
bsize=cdata.blocksize;
rinitial=cdata.randominitial;

%% 初始化
%确定生成多少个基向量（或基向量矩阵）
numbase=ceil(min([NEQ,3*numeig])/bsize);%系数3是为了尽量保证低阶特征对求解的准确性，最后的lanczos矩阵T为numofbase*bsize的方阵
T=zeros(numbase*bsize,numbase*bsize);%先安排空间
Q=zeros(NEQ,numbase*bsize);%QtAQ=T
%生成初始迭代向量
if rinitial==0 %%利用前bsize个单位向量作为初始迭代基底
    Xnow=eye(bsize);
    Xnow(NEQ,bsize)=0;
elseif rinitial==1
    %%这里采用随机子空间迭代法生成初始lanczos迭代向量
    randmatrix=randn(NEQ,bsize);
    Y0=ColSol4eig(randmatrix);
    [Qold,~]=qr(Y0,0);
    for i=1:bsize %这里循环的次数随意
        Ynew=ColSol4eig(Qold);
        [Qnew,~]=qr(Ynew,0);
        Qold=Qnew;
    end
        Xnow=Qold;    
end

invAXnow=ColSol4eig(Xnow);
Mnow=Xnow'*invAXnow;
T(1:bsize,1:bsize)=Mnow;
Q(:,1:bsize)=Xnow;
Xold=zeros(NEQ,bsize);
Bold=zeros(bsize,bsize);

%% 块lanczos变换法的循环过程，bsize=1即为标准的lanczos变换法
for k = 1:numbase-1
    %块状的QR分解
    Rk=invAXnow-Xnow*Mnow-Xold*Bold';
    [Xnew,Bnow]=qr(Rk,0);%%缩减版本的QR分解，见MATLAB程序说明，Xnew的列是单位正交的
    
if reorjudge==1 && bsize==1   
    %这里补上重正交化准则与重正交化过程，对Xnew重正交化
    TK=T(1:k*bsize,1:k*bsize);% k*bsize阶
    beta=norm(Bnow,1);
    normTK=norm(TK,1);
    [V,~]=eig(TK);
    Xnewbackfile=Xnew;
    for num=1:k*bsize
        if beta*abs(V(k*bsize,num))<emsil*normTK
            Y=Q(:,1:k*bsize)*V(:,num);%V的列是单位正交的，因此Y还是单位向量
            Xnew=Xnew-Y*(Y'*Xnew);%选择性重正交化，得出的向量组与Y正交，因此它们张成的空间与Y正交，在QR分解一遍即可
        end
    end
    
    if norm(Xnew-Xnewbackfile)>0 %表明进行了选择性重正交化
        RKnew=Xnew*Bnow;
        [Xnew,Bnow]=qr(RKnew,0);
        %重新对新的Xnew与Bnow张成的RK空间QR分解，在bsize=1时退化到教材P92页的操作
        %块lanczos的求Bnow的做法只是一个猜测，没有理论依据，后续若效果不好可以不更新块lanczos的Bnow
    end
end   

if reorjudge==1 && bsize>1
    %%必须重正交化,而且没调研到什么好的文献，必须全部重正交化
    for numm=1:k
        QQ=Q(:,(numm-1)*bsize+1:numm*bsize);
        Rk=Rk-QQ*(QQ'*Rk);
    end
[Xnew,Bnow]=qr(Rk,0);%%缩减版本的QR分解，见MATLAB程序说明，Xnew的列是单位正交的        
end

    %重正交化完毕
    invAXnew=ColSol4eig(Xnew);
    Mnew=Xnew'*invAXnew;
    
    %填充
    T(k*bsize+1:(k+1)*bsize,k*bsize+1:(k+1)*bsize)=Mnew;
    T(k*bsize+1:(k+1)*bsize,(k-1)*bsize+1:k*bsize)=Bnow;
    T((k-1)*bsize+1:k*bsize,k*bsize+1:(k+1)*bsize)=Bnow';
    Q(:,k*bsize+1:(k+1)*bsize)=Xnew;    
    
    %下一块
    Xold=Xnow;
    Bold=Bnow;
    Xnow=Xnew;
    Mnow=Mnew;
    invAXnow=invAXnew;
    
end

%% 后处理，得到物理空间的正则振型
%T已经是小规模矩阵，直接用MATLAB的eig函数求特征值
[eigT,lamdaT]=eig(T);%%eigT在MATLAB中默认是单位正交向量，这样转换回物理空间自然满足PHIt*M*PHI为单位阵
lamdaT=diag(lamdaT);%成了储存对角元的列向量
%lamdaT要取绝对值（正定阵即为值）最大的numeig个
%y=(linspace(1,numbase*bsize,numbase*bsize))';

lamdasort=sort(lamdaT,'descend');%降序排序，只要前numeig个特征对

lamdaout=zeros(1,numeig);
PHIout=zeros(NEQ,numeig);
mm = M.^(-1/2);
for num=1:numeig
    lamdaout(num)=1/lamdasort(num);%逆迭代法，取倒数
    qq=lamdaT-lamdasort(num);
    [~,hh]=min(abs(qq));
    PHIout(:,num)=mm.*(Q*eigT(:,hh)); %mm.转化为物理空间中的正则振型
end
omegaout=lamdaout.^(1/2);

%%如果有重根，解出的特征值可能有微小的虚部，建议最好开启rinitial以尽可能减少虚部
lamdaout=real(lamdaout);
omegaout=real(omegaout);
PHIout=real(PHIout);

sdata.lamda = lamdaout;%%将处理完的前numeig阶特征对、固有频率返回到全局变量中
sdata.omega = omegaout;%%注意omegaout是圆频率
sdata.PHI = PHIout;

    %%输出各阶主振型
    eigout=fopen('.\Data\EIGEN.OUT', 'w');%%新建一个输出文件
    fprintf(eigout, ['\n' ...
        '      NUMBER OF EIGEN VALUES . . . . . . . . . . (numofeig)  = %10d \n'],...
        numeig);

    ID = sdata.ID;
    NUMNP = cdata.NUMNP;
    for ii=1:numeig
        
    DIS = PHIout(:,ii); 
    fprintf(eigout, '\n\n O M E G A %18.6e', omegaout(ii));
    fprintf(eigout, ['\n\n D I S P L A C E M E N T S\n' ...
        '\n       NODE           X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT\n']);

    D = zeros(3, 1, 'double');
    for II = 1:NUMNP
        D(:) = 0;
        if (ID(1, II) ~= 0) D(1) = DIS(ID(1, II)); end
        if (ID(2, II) ~= 0) D(2) = DIS(ID(2, II)); end
        if (ID(3, II) ~= 0) D(3) = DIS(ID(3, II)); end
        fprintf(eigout, ' %10d        %18.6e%18.6e%18.6e\n', II, D(1), D(2), D(3));
    end

    end

end
