%* *****************************************************************
%* - Function of STAPMAT in solver phase                           *
%*                                                                 *
%* - Purpose:                                                      *
%*     To calculate stresses of P4 element                         *
%*                                                                 *
%* - Call procedures: None                                         *
%*                                                                 *
%* - Called by :                                                   *
%*     SRC/Solver/GetStress.m                                      *
%*                                                                 *
%* *****************************************************************

function P4Stress(NUM, NG)

% Get global data
global cdata;
global sdata;
EC=cdata.EC;
IOUT = cdata.IOUT;
NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ;
E = sdata.E; nu = sdata.nu; t = sdata.t; LM = sdata.LM;
U = sdata.DIS(:, NUM); GC3 = sdata.GC3;

xyz_int = zeros(3,9*NUME); % 9 Gauss Points
strain = zeros(3,9*NUME);
stress = zeros(3,9*NUME);
index = 1;

for N = 1:NUME
    MTYPE = MATP(N); 
    nuN = nu(MTYPE);
    EN = E(MTYPE);
    tN = t(MTYPE);

    DN=EN*tN^3/12/(1-nuN^2)*[1,nuN,0;nuN,1,0;0,0,(1-nuN)/2]; % Constitutive Matrix of plate
    
    % calculate displacements of element nodes
    dis = zeros(12, 1, 'double'); %[w1;tx1;ty1; ... w4;tx4;ty4;]
    for i = 1:12
        if(LM(i,N)>0)
            dis(i) = U(LM(i,N));
        end
    end

    xyz = XYZ(:, N);
    Xs = [xyz(1);xyz(4);xyz(7);xyz(10)];
    Ys = [xyz(2);xyz(5);xyz(8);xyz(11)];
    for i = 1:3
        for j = 1:3
            if EC==0
                Be=C0B(Xs,Ys,GC3(i),GC3(j));
            elseif EC==1
                Be=C1B(Xs,Ys,GC3(i),GC3(j));
            end
            [xx,yy] = RealCoor(Xs,Ys,GC3(i),GC3(j));
            xyz_int(:,index) = [xx; yy; 0.0];
            strain(:,index) = Be * dis; % (3,1) strain on Gauss point
            stress(:,index) = DN * strain(:,index); % (3,1) stress on Gauss point
            index = index + 1;
        end
    end
end

fprintf(IOUT, ['\n\n  S T R A I N  C A L C U L A T I O N S  F O R  ' ...
    'E L E M E N T  G R O U P %4d\n' ... 
    '            POINT  COORDINATE                                 STRAIN\n' ...
    '         X                  Y                  Z                 KX             KY            KXY\n'], NG);

for index=1:9*NUME
   fprintf(IOUT, '%5e  ',xyz_int(:,index)');
   fprintf(IOUT, '   %13.6e',strain(:,index)');
   fprintf(IOUT, '\n'); 
end

fprintf(IOUT, ['\n\n  S T R E S S  C A L C U L A T I O N S  F O R  ' ...
    'E L E M E N T  G R O U P %4d\n\n' ... 
    '            POINT  COORDINATE                                 STRESS\n' ...
    '         X                  Y                  Z                 MX             MY            MXY\n'], NG);

for index = 1:9*NUME
    fprintf(IOUT, '%5e  ',xyz_int(:,index)');
    fprintf(IOUT, '   %13.6e',stress(:,index)');
    fprintf(IOUT, '\n');
end

end

function B=C0B(alpha1,beta1,upsilon,eta)
%     N=[(1+upsilon)*(1+eta),(1-upsilon)*(1+eta),(1-upsilon)*(1-eta),(1+upsilon)*(1-eta)]/4; %C0形函数 w, θ1，θ2独立插值
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
    B=[kx;ky;kxy];
end

function B=C1B(Xs,Ys,kxi,eta)
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
end

function [X,Y]=RealCoor(Xs,Ys,kxi,eta)
    X=kxi*(Xs(1)-Xs(2))/2+(Xs(1)+Xs(2))/2;
    Y=eta*(Ys(1)-Ys(4))/2+(Ys(1)+Ys(4))/2;
end