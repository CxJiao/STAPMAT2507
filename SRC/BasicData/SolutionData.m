%* *****************************************************************
%* - Basic data class of STAPMAT                                   *
%*                                                                 *
%* - Purpose:                                                      *
%*     Storing variables used in solving process                   *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.20                *
%*                                                                 *
%* *****************************************************************

classdef SolutionData
    properties (Constant)
        % Gauss coord, 1D to 3D
        GC1 = double(0.0);
        GC2 = double([sqrt(1/3),-sqrt(1/3)]);
        GC3 = double([sqrt(0.6), 0.0, -sqrt(0.6)]);
        % Gauss weight, 1D to 3D
        GW1 = double(2.0);
        GW2 = double([1.0, 1.0]);
        GW3 = double([5.0/9.0, 8.0/9.0, 5.0/9.0]);
    end
    properties
        % Basic data
        ID;       % int, ID(3, NUMNP), Boundary condition codes (0=free, 1=deleted)
        IDOrigin; % int, backups of ID after computing of NEQ
        X;        % double, X(NUMNP), X coordinates
        Y;        % double, Y(NUMNP), Y coordinates
        Z;        % double, Z(NUMNP), Z coordinates
        R;        % double, R(NEQ, NLCASE), Load vector
        NOD;      % int, NOD(NLOAD), Node number to which this load is applied (1~NUMNP)
        IDIRN;    % int, IDIRN(NLOAD), Degree of freedom number for this load component
                  %                     1 : X-direction;
                  %                     2 : Y-direction;
                  %                     3 : Z-direction;
        FLOAD;    % double, FLOAD(NLOAD), Magnitude of load
        PRSID;    % int, PRSID(NUME), ID of pressure load applied on each plate elements
        PRESS;    % double, PRESS(NSTEPS,max(PRSID)), time-variant pressure load of different PRSID
        
        % Element data
        NUME;     % int, number of elements for TRUSS & P4
        NUME2;    % int, number of elements for H20
        NNODE;    % int, number of nodes in an element
        NINIP;    % int, number of integration points in an element
        NDOF;     % int, the DOF of displacement
        NSTIFF;   % int, the number of number in element stiffness matrix
        XYZ;      % double, XYZ(3*NNODE, NUME), element position for TRUSS &P4
        XYZ2;     % double, XYZ(3*NNODE, NUME2), element position for H20
        ELNOD     % int, ELNOD(NNODE,NUME)
        InitCoord;  % double array, integration coordinates
        InitWeight; % double array, integration weights
        
        % Material data
        NUMMAT;     % int, the number of types of material 
        E;          % double array, Young's Modulus for TRUSS & P4
        E2;         % double array, Young's Modulus for H20
        nu;         % double array, possion ratio
        AREA;       % double array, cross-sectional constants
        t;          % double array, plate thickness for P4 element
        rho;        % double array, density
        MATP;       % int, MATP(NUME), types of elements for TRUSS & P4
        MATP2;      % int, MATP(NUME), types of elements for H20
        
        % Solve data
        NEQ;      % int, Number of equations
        NWK;      % Number of matrix elements
        MK;       % Maximum half bandwidth
        MHT;      % int, MHT(NEQ), Vector of column heights
        LM;       % int, LM(6, NUME), Connectivity matrix for TRUSS
        LM2;      % int, LM(60, NUME), Connectivity matrix for H20
        MAXA;     % int, MAXA(NEQ+1)，K的各对角元在A中的地址向量，一共有NEQ+1个元素，第一个元素为1，第二个元素为2
        STIFF;    % double, STIFF(NWK), store the elements of stiffness matrix
        MASS;     % double, MASS(NEQ), store the diag elements of lumped mass matrix，是一个列阵
        MASSC;    % double, MASSC(NWK), store the elements of consistent matrix
        stiff4eig;%用来算广义特征对的刚度阵，生成时与STIFF同时生成，但进行M^(-1/2）*stiff4eig*M^(-1/2）后再LDLT分解，使问题变为一个标准特征值问题
        lamda;    %用来储存广义特征值的列向量，从小到大储存。lamda(1,numofeig)
        omega;    %用来储存固有频率的列向量，为sqrt(lamda)。omega(1,numofeig)，从左至右与与lamda对应
        PHI;      %用来储存正则振型的矩阵，为PHI(NEQ，numofeig)，从左至右与lamda对应

        % Dynamic solve data
        ENDTIME;  % double
        CparaM;   % double, Rayleigh Damping Matrix - parameter with Mass Matrix
        CparaK;   % double, Rayleigh Damping Matrix - parameter with Stiffness Matrix
        DyAlphaF; % double, parameters of generalized alpha method
        DyAlphaM;
        DyGamma;
        DyBeta;
        TimeStep; % double
        MassType; % int, type of mass matrix(1-lumped,2-consistent) in dynamics
        NSTEPS;   % int, amounts of Time Steps <- fixed TimeStep
        Ddis;     % double(NEQ), initial and updated displacements in dynamic solutions
        Dvol;     % double(NEQ), initial and updated volecities in dynamic solutions
        Dacc;     % double(NEQ), initial and updated accelerations in dynamic solutions
        EnKine;   % double(NSTEPS), kinetic energy
        InitialDis;% int, initial displacements type
        NVNL;     % int, amount of Time-variant Nodal loads
        DLDC;     % int(NVNL), Load case of Time-variant Nodal load
        DNOD;     % int(NVNL), Node number to which this load is applied (1~NUMNP)
        DDIRE;    % int(NVNL), Degree of freedom number for this load component
                  %                     1 : X-direction;
                  %                     2 : Y-direction;
                  %                     3 : Z-direction;
        DVNL;     % double(NSTEPS,NVNL), magtitude of Time-variant Load

        % Result data
        DIS;      % double, DIS(NEQ, NLCASE), Displacement of nodes
        STRAIN;   % double, STRAIN(NEQ, NLCASE), Strain
        STRESS;   % double, STRESS(NEQ, NLCASE), Stress
        TDIS;     % cell, TDIS{NLCASE*[NEQ,NSTEP]}, Time-variant Displacement of node
        TDISW;    % cell, TDIS{NLCASE*[NUMNP,NSTEP]}, Time-variant Displacement w of node (for P4)
        TSTRN;     % cell, TDIS{NLCASE*[NUME,NSTEP]}, Time-variant average stress of element
        TSTRS;     % cell, TDIS{NLCASE*[NUME,NSTEP]}, Time-variant average strain of element
    end
end