%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     Compute the stiffness matrix of P4 element                  *
%*     For Kirchhoff-Love thin plate                               *
%*                                                                 *
%* - Call procedures:                                              *
%*                                                                 *
%* - Called by :                                                   *
%*                                                                 *
%* *****************************************************************

function P4Stiff(NG)

% Init variables of the element
InitP4();

% Read Material and Elements
ReadP4(NG);

fprintf('Solution phase ...\n\n');

% calculate addresses of diagonal elements
Addres(NG);

% Data check Or Solve
global cdata;
cdata.EC = 1; % 0=C0stiff 1=C1stiff. default is C1
if (cdata.MODEX == 0) 
    cdata.TIM(3,:) = clock;
    cdata.TIM(4,:) = clock;
    cdata.TIM(5,:) = clock;
    return; 
end

% Assemble structure stiffness matrix
Assemble(NG);

end

function InitP4()
global sdata;
sdata.NNODE = 4; % 4 nodes per element
sdata.NDOF = 3;  % 3 degree of freedom per node. w, thetax, thetay
end

function Assemble(NG)
global sdata;
global cdata;
% Initial STIFF MASSC
if( NG == 1 )
    sdata.STIFF = zeros(sdata.NWK, 1, 'double');
    if(sdata.MassType == 2)
        sdata.MASSC = zeros(sdata.NWK, 1, 'double');
    end
end
EC=cdata.EC;
NUME = sdata.NUME; MATP = sdata.MATP;
E = sdata.E;nu = sdata.nu; t = sdata.t; rho = sdata.rho;
XYZ = sdata.XYZ; LM = sdata.LM;
GC3 = sdata.GC3; GW3 = sdata.GW3;
GC1 = sdata.GC1; GW1 = sdata.GW1;

for N = 1:NUME
% Stiff Matrix
    ES=zeros(12,12, 'double'); % Stiff matrix of one element
    MTYPE = MATP(N); 
    nuN = nu(MTYPE); EN = E(MTYPE); tN = t(MTYPE); rhoN=rho(MTYPE);
    DN=EN*tN^3/12/(1-nuN^2)*[1,nuN,0;nuN,1,0;0,0,(1-nuN)/2]; % Constitutive Matrix
    XN=[XYZ(1,N);XYZ(4,N);XYZ(7,N);XYZ(10,N)];
    YN=[XYZ(2,N);XYZ(5,N);XYZ(8,N);XYZ(11,N)];
    
    if EC==0
         %对Kb采用两点高斯积分，对Ks采用单点高斯积分
        Kb=zeros(12,12, 'double');
        for i = 1:3
            for j = 1:3
                [Je,BDB,~]=C0stiff(DN,XN,YN,GC3(i),GC3(j),1000*EN);
                Je=det(Je);
                Kb=Kb+GW3(i)*GW3(j)*Je*BDB;
            end
        end
        [Je,~,BDBs]=C0stiff(DN,XN,YN,GC1(1),GC1(1),100*EN);
        Je=det(Je);
        Ks=GW1(1)*GW1(1)*Je*BDBs;
        ES=Ks+Kb;
    elseif EC==1
        for i = 1:3
            for j = 1:3
                [Je,BDB]=C1stiff(DN,XN,YN,GC3(i),GC3(j));
                Je=det(Je);
                ES=ES+GW3(i)*GW3(j)*Je*BDB;
            end
        end
    end
    %   SRC/Mechanics/ADDBAN.m Assemble the structure stiff matrix
    ADDBAN(ES, LM(:, N));
    
% Mass Matrix
    Aera=abs(XN(1)-XN(2))*abs(YN(1)-YN(4));
    WN = rhoN * tN * Aera; % total mass of element
    EM=(WN/4)*[1;0;0;1;0;0;1;0;0;1;0;0];
    ADDMASS(EM, LM(:, N));
    
    if(sdata.MassType == 2)
        EMC=zeros(12,12, 'double');
        if EC==0
            for i = 1:3
                for j = 1:3
                    [Je,Ne]=C0shape(XN,YN,GC3(i),GC3(j));
                    NTN=Ne'*Ne;
                    Je=det(Je);
                    EMC=EMC+GW3(i)*GW3(j)*Je*NTN;
                end
            end
        elseif EC==1
            for i = 1:3
                for j = 1:3
                    [Je,Ne]=C1shape(XN,YN,GC3(i),GC3(j));
                    NTN=Ne'*Ne;
                    Je=det(Je);
                    EMC=EMC+GW3(i)*GW3(j)*Je*NTN;
                end
            end
        end
        EMC= rhoN * tN * EMC;
        ADDMASSC(EMC, LM(:, N));
    end
    
% Pressure Assemble
    PE=zeros(12,1,'double');
    if EC==0
        for i = 1:3
            for j = 1:3
                [Je,Ne]=C0shape(XN,YN,GC3(i),GC3(j));
                Je=det(Je);
                PE=PE+GW3(i)*GW3(j)*Je*Ne(1,:)';
            end
        end
    elseif EC==1
        for i = 1:3
            for j = 1:3
                [Je,Ne]=C1shape(XN,YN,GC3(i),GC3(j));
                Je=det(Je);
                PE=PE+GW3(i)*GW3(j)*Je*Ne';
            end
        end
    end
    ADDPRESS(PE, N);
end

% The third time stamp
cdata.TIM(3, :) = clock;

end

function [Je,BDB,BDBs]=C0stiff(D,alpha1,beta1,upsilon,eta,G)
    % Constructed from Mindlin Plate element
    % alpha1=[x1;x2;x3;x4], beta1=[y1;y2;y3;y4]. Absolute location
    % the first point is (1,1) on the up-right corner
    % upsilon=x of Gauss point, eta=y of Gauss point in the range of [-1,1]

    N=[(1+upsilon)*(1+eta),(1-upsilon)*(1+eta),(1-upsilon)*(1-eta),(1+upsilon)*(1-eta)]/4; %C0形函数 w, θ1，θ2独立插值
    pNpeta=[(1+eta),-(1+eta),-(1-eta),(1-eta); % dN/dupsilon
        (1+upsilon),(1-upsilon),-(1-upsilon),-(1+upsilon)]/4; % dN/deta
    Je=pNpeta*[alpha1,beta1]; %2x2 [dX/dupsilon,dY/dupsilon;dX/deta,dY/deta]
    Je1=inv(Je);
    pNpalpha=Je1*pNpeta;
    N1alpha=pNpalpha(1,1)*diag([1,1,1]);
    N2alpha=pNpalpha(1,2)*diag([1,1,1]);
    N3alpha=pNpalpha(1,3)*diag([1,1,1]);
    N4alpha=pNpalpha(1,4)*diag([1,1,1]);
    N1beta=pNpalpha(2,1)*diag([1,1,1]);
    N2beta=pNpalpha(2,2)*diag([1,1,1]);
    N3beta=pNpalpha(2,3)*diag([1,1,1]);
    N4beta=pNpalpha(2,4)*diag([1,1,1]);
    Z0=zeros(3,3);
    I=diag(ones(3,1));
    dr1=[N1alpha,N2alpha,N3alpha,N4alpha]; % 仅代表系数矩阵，需乘结点位移才是真正的[dw/dx;dthetax/dx;dthetay/dx]
    dr2=[N1beta,N2beta,N3beta,N4beta]; % 仅代表系数矩阵，需乘结点位移才是真正的[dw/dy;dthetax/dy;dthetay/dy]
    kx=-dr1(2,:);
    ky=-dr2(3,:);
    kxy=-(dr1(3,:)+dr2(2,:));
    %广义应变矩阵dgammar,dkappal,dkappam,detall,detam,dwl,dwm,dlamma,dmu
    Bb=[kx;ky;kxy];
    %刚度矩阵
    BDB=Bb'*D*Bb;
    %罚函数矩阵
    Bs1=[pNpalpha(1,1),-N(1),0
        pNpalpha(2,1),0,-N(1)];
    Bs2=[pNpalpha(1,2),-N(2),0;
        pNpalpha(2,2),0,-N(2)];
    Bs3=[pNpalpha(1,3),-N(3),0;
        pNpalpha(2,3),0,-N(3)];
    Bs4=[pNpalpha(1,4),-N(4),0;
        pNpalpha(2,4),0,-N(4)];
    Bs=[Bs1,Bs2,Bs3,Bs4];
    BDBs=Bs'*G*Bs;
end

function [Je,BDB]=C1stiff(D,Xs,Ys,kxi,eta)
    % Constructed from Kirchhoff Plate element
    % Xs=[x1;x2;x3;x4], Ys=[y1;y2;y3;y4]. Absolute location
    % the first point is (1,1) on the up-right corner
    % (kxi,eta) is the corrnidate of Gauss point in the range of [-1,1]
    a=abs(Xs(1)-Xs(2))/2;
    b=abs(Ys(1)-Ys(4))/2;
    Je=[a,0;0,b];
    B=zeros(3,12);
    kxis=[1,-1,-1,1];  etas=[1,1,-1,-1];                                                                                                                              
    for n=1:4
        B(:,3*n-2:3*n)=[(1+eta*etas(n))*3*kxi*kxis(n)/(4*a^2),0,(1+eta*etas(n))*kxis(n)*(1+3*kxi*kxis(n))/(4*a); ...
            (1+kxi*kxis(n))*3*eta*etas(n)/(4*b^2),-(1+kxi*kxis(n))*etas(n)*(1+3*eta*etas(n))/(4*b),0; ...
            kxis(n)*etas(n)*(-4+3*eta^2+3*kxi^2)/(4*a*b),-(3*eta^2+2*eta*etas(n)-1)*kxis(n)/(4*a),(3*kxi^2+2*kxi*kxis(n)-1)*etas(n)/(4*b)];
    end
    BDB=B'*D*B;
end

function [Je,N]=C0shape(Xs,Ys,upsilon,eta)
    a=abs(Xs(1)-Xs(2))/2;
    b=abs(Ys(1)-Ys(4))/2;
    Je=[a,0;0,b];
    N3=[(1+upsilon)*(1+eta)*diag([1,1,1]),(1-upsilon)*(1+eta)*diag([1,1,1]), ...
        (1-upsilon)*(1-eta)*diag([1,1,1]),(1+upsilon)*(1-eta)*diag([1,1,1])]/4; %C0形函数 w, θ1，θ2独立插值
    N=N3;
end

function [Je,N]=C1shape(Xs,Ys,kxi,eta)
    N=zeros(1,12);
    a=abs(Xs(1)-Xs(2))/2;
    b=abs(Ys(1)-Ys(4))/2;
    Je=[a,0;0,b];
    kxis=[1,-1,-1,1];  etas=[1,1,-1,-1];   
    for n=1:4
        N(3*n-2:3*n)=[(kxis(n)*kxi+1)*(etas(n)*eta+1)*(2+kxis(n)*kxi+etas(n)*eta-kxi^2-eta^2), ...
            b*etas(n)*(kxis(n)*kxi+1)*(etas(n)*eta+1)^2*(etas(n)*eta-1), ...
            -a*kxis(n)*(kxis(n)*kxi+1)^2*(kxis(n)*kxi-1)*(etas(n)*eta+1)];
    end
    N=0.125*N;
%     N=0.125*[(kxi+1)*(eta+1)*(2+eta-eta^2+kxi-kxi^2),b*(eta-1)*(eta+1)^2*(kxi+1),-a*(eta+1)*(kxi-1)*(kxi+1)^2, ...
%         (kxi-1)*(eta+1)*(-2-eta+eta^2+kxi+kxi^2),-b*(eta-1)*(eta+1)^2*(kxi-1),-a*(eta+1)*(kxi-1)^2*(kxi+1), ...
%         -(kxi-1)*(eta-1)*(-2+eta+eta^2+kxi+kxi^2),-b*(eta-1)^2*(eta+1)*(kxi-1),a*(eta-1)*(kxi-1)^2*(kxi+1), ...
%         (kxi+1)*(eta-1)*(-2+eta+eta^2-kxi+kxi^2),b*(eta-1)^2*(eta+1)*(kxi+1),a*(eta-1)*(kxi-1)*(kxi+1)^2];
%     dNdkxi=0.125*[(eta+1)*(3+eta-eta^2-3*kxi^2),b*(eta-1)*(eta+1)^2,-a*(eta+1)*(-1+2*kxi+3*kxi^2), ...
%         (eta+1)*(-3-eta+eta^2+3*kxi^2),-b*(eta-1)*(eta+1)^2,-a*(eta+1)*(-1-2*kxi+3*kxi^2), ...
%         -(eta-1)*(-3+eta+eta^2+3*kxi^2),-b*(eta-1)^2*(eta+1),a*(eta-1)*(-1-2*kxi+3*kxi^2), ...
%         (eta-1)*(-3+eta+eta^2+3*kxi^2),b*(eta-1)^2*(eta+1),a*(eta-1)*(-1+2*kxi+3*kxi^2)];
%     dNdeta=0.125*[(kxi+1)*(3-3*eta^2+kxi-kxi^2),b*(kxi+1)*(-1+2*eta+3*eta^2),-a*(kxi-1)*(kxi+1)^2, ...
%         (kxi-1)*(-3+3*eta^2+kxi+kxi^2),-b*(kxi-1)*(-1+2*eta+3*eta^2),-a*(kxi+1)*(kxi-1)^2, ...
%         -(kxi-1)*(-3+3*eta^2+kxi+kxi^2),-b*(kxi-1)*(-1-2*eta+3*eta^2),a*(kxi+1)*(kxi-1)^2, ...
%         (kxi+1)*(-3+3*eta^2-kxi+kxi^2),b*(kxi+1)*(-1-2*eta+3*eta^2),a*(kxi-1)*(kxi+1)^2];
end

