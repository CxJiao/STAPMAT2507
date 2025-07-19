clc
clear

% Generate the nodes and elemets information of plate mesh
% Node Sequence is Right-up, Left-up, Left-down, Right-down
La=1.6; % Length in x-direction (m)
Lb=1;  % Length in y-direction (m)
Na=20; % element number in La
Nb=20; % element number in Lb
% Boundary Condition type on the 4 edges [AB BC CD DA] 0=free, 1=fix
BCw=[1,1,1,1]; % BC about w
BCtx=[1,0,1,0]; % BC about thetax
BCty=[0,1,0,1]; % BC about thetay
LDP=[0.0,0.0,0.0]; % Coordinate of the Left-down point of the plate (1st Node)
NUME=Na*Nb;
NUMNP=(Na+1)*(Nb+1);
dx=La/Na;
dy=Lb/Nb;

MTYPE=ones(NUME,1);
PRSID=ones(NUME,1);
XYZ=zeros(NUMNP,3);
DOF=zeros(NUMNP,3);
NID=zeros(Na+1,Nb+1);
EID=[];
CP=LDP;
index=1;
for j=1:Nb+1
    for i=1:Na+1
        XYZ(index,:)=CP;
        NID(j,i)=index;
        % manage the boundary condition
        if i == 1
            if j == 1
                NOF=int64([BCw(3)||BCw(2),BCtx(3)||BCtx(2),BCty(3)||BCty(2)]);
            elseif j==Nb+1
                NOF=int64([BCw(1)||BCw(2),BCtx(1)||BCtx(2),BCty(1)||BCty(2)]); 
            else
                NOF=[BCw(2),BCtx(2),BCty(2)];
            end
        elseif i == Na+1
            if j == 1
                NOF=int64([BCw(3)||BCw(4),BCtx(3)||BCtx(4),BCty(3)||BCty(4)]);
            elseif j==Nb+1
                NOF=int64([BCw(1)||BCw(4),BCtx(1)||BCtx(4),BCty(1)||BCty(4)]); 
            else
                NOF=[BCw(4),BCtx(4),BCty(4)];
            end
        else
            if j == 1
                NOF=[BCw(3),BCtx(3),BCty(3)];
            elseif j == Nb+1
                NOF=[BCw(1),BCtx(1),BCty(1)];
            else
                NOF=[0,0,0];
            end
        end
        DOF(index,:)=NOF;
        % manage the element information
        if i > 1 && j > 1
            if i==2 && j==2
                EID=[NID(j,i),NID(j,i-1),NID(j-1,i-1),NID(j-1,i)];
            else
                EID=[EID;[NID(j,i),NID(j,i-1),NID(j-1,i-1),NID(j-1,i)]];
            end
        end
        CP(1)=CP(1)+dx;
        index=index+1;
    end
    CP(1)=0.0;
    CP(2)=CP(2)+dy;
end

TN=(1:1:NUMNP)';
TE=(1:1:NUME)';
NOD_INFO=[TN,DOF,XYZ];
ELE_INFO=[TE,EID,MTYPE,PRSID];
writematrix('**************NODES***************','MESH_INFO.txt');
writematrix(NOD_INFO,'MESH_INFO.txt','Delimiter',' ','WriteMode','append');
writematrix('*************ELEMENTS*************','MESH_INFO.txt','WriteMode','append');
writematrix(ELE_INFO,'MESH_INFO.txt','Delimiter',' ','WriteMode','append');

