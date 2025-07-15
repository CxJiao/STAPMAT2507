% ******************************************************************
%*                                                                 *
%*                         S T A P M A T                           *
%*                                                                 *
%*      An In-CORE SOLUTION STATIC ANALYSIS PROGRAM IN MATLAB      *
%*      Adapted from STAP90 (FORTRAN 90) for teaching purpose      *
%*                                                                 *
%*  Computational Dynamics Group, School of Aerospace Engineering  *
%*  Tsinghua University, 2019.02.20                                *
%*                                                                 *
%* *****************************************************************

% Set paths of functions 
AddPath();

% Define Global Variables   
global cdata;
global sdata;
cdata = ControlData;
sdata = SolutionData;

% Read InPut file

% fname = 'H20.in';
% fname = 'beam.in';
fname = 'Dchair.in';

ReadFile(fname);

% Write basic data of program 
WriteParasOut();

% Form the stiffness matrix
GetStiff();

% 进行特征值的计算
%if cdata.numofeig>=1  %%若numofeig为0，跳过该步
%processstiff4eig();%将问题转化为标准特征值问题
%blocklanczos();%进行块lanczos迭代
%end

% Triangularize stiffness matrix
Solve();

% 进行特征值的计算，静力学之后已经解出了K的LDLT分解，特征值跟静力学一起做
if cdata.numofeig>=1 && cdata.MODEX == 1 %%若numofeig为0，跳过该步，注意一定为静力学模式才能做
processstiff4eig2();%将问题转化为标准特征值问题
blocklanczos();%进行块lanczos迭代
end

% Finalize
Finalize();

% ----------------------- Functions -----------------------------------

% Functions
% Add paths of functions
function AddPath()
clear;
close all;
clc;

% addpath .\SRC\Initiation
% addpath .\SRC\BasicData
% addpath .\SRC\Mechanics
% addpath .\SRC\Mechanics\Truss
% addpath .\SRC\Mechanics\H20
% addpath .\SRC\Solver
addpath ./SRC/Initiation
addpath ./SRC/BasicData
addpath ./SRC/Mechanics
addpath ./SRC/Mechanics/Truss
addpath ./SRC/Mechanics/H20
addpath ./SRC/Solver
end

function Finalize()
global cdata;
TIM = cdata.TIM;
time = zeros(5, 1, 'double');
time(1) = etime(TIM(2,:), TIM(1,:));
time(2) = etime(TIM(3,:), TIM(2,:));
time(3) = etime(TIM(4,:), TIM(3,:));
time(4) = etime(TIM(5,:), TIM(4,:));
time(5) = etime(TIM(5,:), TIM(1,:));
% etime衡量消耗时间,新版本推荐改写为etime
% 这里最好在特征值算完之后再填进去一个时刻，好看看特征值计算单个的效率，补进去TIM(4,:)，之前的后移一位

fprintf(cdata.IOUT, ['\n\n' ...
    ' S O L U T I O N   T I M E   L O G   I N   S E C\n\n' ...
    '     TIME FOR INPUT PHASE  . . . . . . . . . . . . . . = %12.2f\n' ...
    '     TIME FOR CALCULATION OF STIFFNESS MATRIX  . . . . = %12.2f\n' ...
    '     TIME FOR FACTORIZATION OF STIFFNESS MATRIX  . . . = %12.2f\n' ...
    '     TIME FOR LOAD CASE SOLUTIONS  . . . . . . . . . . = %12.2f\n\n' ...
    '      T O T A L   S O L U T I O N   T I M E  . . . . . = %12.2f\n'], ...
    time(1), time(2), time(3), time(4),time(5));

fprintf(['\n' ...
    ' S O L U T I O N   T I M E   L O G   I N   S E C\n\n' ...
    '     TIME FOR INPUT PHASE  . . . . . . . . . . . . . . = %12.2f\n' ...
    '     TIME FOR CALCULATION OF STIFFNESS MATRIX  . . . . = %12.2f\n' ...
    '     TIME FOR FACTORIZATION OF STIFFNESS MATRIX  . . . = %12.2f\n' ...
    '     TIME FOR LOAD CASE SOLUTIONS  . . . . . . . . . . = %12.2f\n\n' ...
    '      T O T A L   S O L U T I O N   T I M E  . . . . . = %12.2f\n'], ...
    time(1), time(2), time(3), time(4),time(5));

fclose(cdata.IIN);
fclose(cdata.IOUT);
end