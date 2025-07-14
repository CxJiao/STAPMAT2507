%* *****************************************************************
%* - Basic data class of STAPMAT                                   *
%*                                                                 *
%* - Purpose:                                                      *
%*     Storing variables which control the running of STAPMAT      *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.20                *
%*                                                                 *
%* *****************************************************************
classdef ControlData
    properties
        NUMNP;         % Total number of nodal points
                       % = 0 : Program stop

        NPAR;          %   int(3, NUMEG), Element group control data
                       %   NPAR(1, :) - Element type
                       %             1 : Truss element
                       %             2 : H20 element
                       %   NPAR(2, :) - Number of elements
                       %   NPAR(3, :) - Number of different sets of material
                       %             and cross-sectional constants
        NUMEG;         % Total number of element groups, > 0
        NLCASE;        % Number of load case (>0)
        LL;            % Load case number
        NLOAD;         % The number of concentrated loads applied in this load case

        MODEX;         % Solution mode: 0 - data check only; 1 - execution

        TIM;           % Timing information
        HED;           % Master heading information for usr in labeling the output

        IIN;           % file pointer used for input
        IOUT;          % file pointer used for output
        
        numofeig;      %所需的特征对个数（小于总方程数，若不想求特征对，这里输入0）
        blocksize;     %运用块lanczos迭代的块尺寸
        randominitial; %是否用随机方法生成迭代初始向量（1是0否）
    end
end
