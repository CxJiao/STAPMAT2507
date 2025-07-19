%* *****************************************************************
%* - PostProcess of STAPMAT after Solving                          *
%*                                                                 *
%* - Called by :                                                   *
%*     stapmat.m                                                   *
%*                                                                 *
%* *****************************************************************

function PostProcess(FILE,DATATYPE)
global cdata
global sdata
FNAME=strcat('./Data/',FILE);
cdata.IPOST=fopen(FNAME, 'w');
TITLE=cdata.HED; 
fprintf(cdata.IPOST, 'TITLE = "%s"\n', TITLE);
sdata.TDISW=cell(cdata.NLCASE);
sdata.TSTRN=cell(cdata.NLCASE);
sdata.TSTRS=cell(cdata.NLCASE);
for L = 1:cdata.NLCASE
    for N = 1:cdata.NUMEG
        NPAR1 = cdata.NPAR(1,N);
        if (NPAR1 == 3) 
            P4Conclude(L);
            P4Post(L,DATATYPE);
        end
    end
end

end

function P4Conclude(NLC)
global sdata
global cdata
EC=cdata.EC;
NUMNP = cdata.NUMNP; NUME = sdata.NUME; NSTEPS=sdata.NSTEPS;
E = sdata.E; nu = sdata.nu; t = sdata.t; LM = sdata.LM;
XYZ=sdata.XYZ;ID=sdata.ID; MATP = sdata.MATP; GC3 = sdata.GC3;
TDIS=sdata.TDIS{NLC};

TDISW=zeros(NUMNP,NSTEPS);
TSTRN=zeros(3*NUME,NSTEPS); % [kx;ky;kxy;...]
TSTRS=zeros(3*NUME,NSTEPS); % [mx;my;mxy;...]

for I=1:NSTEPS
    for J=1:NUMNP
        % manage the displacement w
        EQ=ID(1,J); 
        if EQ >0
            TDISW(J,I)=TDIS(EQ,I);
        else
            TDISW(J,I)=0;
        end
    end
    
    for N = 1:NUME
        MTYPE = MATP(N);
        nuN = nu(MTYPE);
        EN = E(MTYPE);
        tN = t(MTYPE);
        
        DN=EN*tN^3/12/(1-nuN^2)*[1,nuN,0;nuN,1,0;0,0,(1-nuN)/2]; % Constitutive Matrix of plate
        
        % calculate displacements of element nodes
        dis = zeros(12, 1, 'double'); %[w1;tx1;ty1; ... w4;tx4;ty4;]
        Strain=zeros(3,1,'double');
        Stress=zeros(3,1,'double');
        for i = 1:12
            EQ=LM(i,N);
            if EQ>0
                dis(i) = TDIS(EQ,I);
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
                strain = Be * dis; % (3x1) strain on Gauss points
                stress = DN * strain; % (3x1) stress on Gauss points
                Strain = Strain + (1/9).*strain; % average of 9 points
                Stress = Strain + (1/9).*stress; % average of 9 points
            end
        end
        TSTRN(3*N-2:3*N,I)=Strain;
        TSTRS(3*N-2:3*N,I)=Stress;
    end
    
end
sdata.TDISW{NLC}=TDISW;
sdata.TSTRN{NLC}=TSTRN;
sdata.TSTRS{NLC}=TSTRS;
end

function P4Post(NLC,DATATYPE)
global sdata
global cdata

IPOST=cdata.IPOST;
NUMNP = cdata.NUMNP; NUME = sdata.NUME; 
ELNOD=sdata.ELNOD; NSTEPS=sdata.NSTEPS;
TDISW=sdata.TDISW{NLC}; % displacement result
X=sdata.X; Y=sdata.Y; Z=sdata.Z;
dt=sdata.TimeStep;
time=0;
for K=1:length(DATATYPE)
    if DATATYPE(K) == 'U'
        fprintf(IPOST, 'VARIABLES = "X", "Y", "Z", "Displacement_X", "Displacement_Y", "Displacement_Z"\n\n');
        for I=1:NSTEPS
            if I==1
                fprintf(IPOST, 'Zone T="Time= %.3e",F=FEPOINT, N= %d, E= %d, ET=QUADRILATERAL\n' ...
                    ,time, NUMNP, NUME);
            else
                fprintf(IPOST, 'Zone T="Time= %.3e",F=FEPOINT, N= %d, E= %d, ET=QUADRILATERAL, D=(FECONNECT)\n' ...
                    ,time, NUMNP, NUME);
            end
           for J=1:NUMNP
               fprintf(IPOST, '%.3f  %.3f  %.3f  %.3e  %.3e  %.3e \n', ...
                   X(J),Y(J),Z(J),0.0,0.0,TDISW(J,I));
           end
           if I==1
               for P=1:NUME
                   fprintf(IPOST,'%d  %d  %d  %d\n',ELNOD(1,P),ELNOD(2,P),ELNOD(3,P),ELNOD(4,P));
               end
           end
           fprintf(IPOST, '\n');
           time = time+dt;
        end
    end
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
