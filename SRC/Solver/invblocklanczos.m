function [T,Q] = invblocklanczos(A,m,r)
%BLOCKLANCZOS QTAQ=T，其中T中的小矩阵为m阶的，要迭代r次，即最后的T矩阵是m*r阶的方阵，逆迭代，只是个示例程序
%没有被总程序调用，采用了随机子空间迭代生成初始lanczos迭代阵，小矩阵测试效果很好，T的阶数m*r大约取3p
%p是需要的特征对数
invA=inv(A);

n=size(A,1);
T=zeros(m*r,m*r);
Q=zeros(n,m*r);
%X1=eye(m);%%初始化
%X1(n,m)=0;
omega=randn(n,m);
Y0=invA*omega;
[Qold,~]=qr(Y0,0);
for i=1:m-1
    Ynew=invA*Qold;
    [Qnew,~]=qr(Ynew,0);
    Qold=Qnew;
end
    X1=Qold;

M1=X1'*invA*X1;
T(1:m,1:m)=M1;
Q(:,1:m)=X1;

Xold=zeros(n,m);
Bold=zeros(m,m);
Xnow=X1;
Mnow=M1;
    for k=1:r-1
        %%块状的QR分解
        Rk=invA*Xnow-Xnow*Mnow-Xold*Bold';
        [Xnew,Bnow]=qr(Rk);
        Xnew=Xnew(:,1:m);
        Bnow=Bnow(1:m,1:m);
        Mnew=Xnew'*invA*Xnew;

        %%填充
        T(k*m+1:(k+1)*m,k*m+1:(k+1)*m)=Mnew;
        T(k*m+1:(k+1)*m,(k-1)*m+1:k*m)=Bnow;
        T((k-1)*m+1:k*m,k*m+1:(k+1)*m)=Bnow';
        Q(:,k*m+1:(k+1)*m)=Xnew;
        
        %%下一步
        Xold=Xnow;
        Bold=Bnow;
        Xnow=Xnew;
        Mnow=Mnew;
    end
end

