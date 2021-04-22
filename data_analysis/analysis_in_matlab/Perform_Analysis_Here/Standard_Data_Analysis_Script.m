%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled  
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  nick.battista@unc.edu
% Date Created: May 27th, 2015
% Institution: UNC-CH
%
% This code is capable of creating Lagrangian Structures using:
% 	1. Springs
% 	2. Beams (*torsional springs)
% 	3. Target Points
%	4. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-Tension)")
%
% One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting lengths, etc
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
% 
% If you would like us to add a specific muscle model, please let Nick (nick.battista@unc.edu) know.
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: This is where you can perform data analysis of a simulation
%           from your stored viz_IB2d information
%
%     Note: 
%           (1) USER-DEFINED functions should be made made to analyze specific 
%               data sets, regions of interest, etc. (see
%               EXAMPLE_FOR_DATA_ANALYSIS for an example of this)
%           (2) MUST make sure to 'addpath' to where DA_Blackbox is, i.e.,
%               line 64
%           (3) MUST make sure to set path to desired dataset, i.e., in line
%               59 and 60
%           (4) EASIEST to move your viz_IB2d and hier_IB2d_data folders
%               into THIS folder
%          
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Standard_Data_Analysis_Script()

% TEMPORAL INFO FROM input2d %
dt = 2.5e-5;      % Time-step
Tfinal = 0.05;   % Final time in simulation
pDump=50;       % Note: 'print_dump' should match from input2d

% DATA ANALYSIS INFO %
start=25;                             % 1ST interval # included in data analysis
finish=25;                           % LAST interval # included in data analysis 
dump_Times = (start:1:finish)*pDump; % Time vector when data was printed in analysis

% SET PATH TO DESIRED viz_IB2d / hier_IB2d_data DATA %
pathViz = 'viz_IB2d';
pathForce= 'hier_IB2d_data';

% SET PATH TO DA_BLACKBOX %
addpath('../DA_Blackbox');  % <--- MAKE SURE THIS IS PASSED ACCORDINGLY


for i=start:1:finish
        
        % Points to desired data viz_IB2d data file
        if i<10
           numSim = ['000', num2str(i) ];
        elseif i<100
           numSim = ['00', num2str(i) ];
        elseif i<1000
           numSim = ['0', num2str(i)];
        else
           numSim = num2str(i);
        end

        % Imports immersed boundary positions %
        [xLag,yLag] = give_Lag_Positions(pathViz,numSim);

        % Imports (x,y) grid values and ALL Eulerian Data %
        %                      DEFINITIONS 
        %          x: x-grid                y: y-grid
        %       Omega: vorticity           P: pressure
        %    uMag: mag. of velocity  
        %    uX: mag. of x-Velocity   uY: mag. of y-Velocity  
        %    U: x-directed velocity   V: y-directed velocity
        %    Fx: x-directed Force     Fy: y-directed Force
        %
        %  Note: U(j,i): j-corresponds to y-index, i to the x-index
        %
        % 
        % <<- CHOOSE WHAT EULERIAN DATA YOU WANT TO ANALYZE ->>
        %
        Eulerian_Flags(1) = 0;   % OMEGA
        Eulerian_Flags(2) = 0;   % PRESSURE
        Eulerian_Flags(3) = 0;   % uMAG
        Eulerian_Flags(4) = 0;   % uX (mag. x-component of velocity)
        Eulerian_Flags(5) = 0;   % uY (mag. x-component of velocity)
        Eulerian_Flags(6) = 0;   % uVEC (vector components of velocity: U,V)
        Eulerian_Flags(7) = 0;   % Fx (x-component of force )
        Eulerian_Flags(8) = 0;   % Fy (y-component of force)
        Eulerian_Flags(9) = 0;   % C (concentration)
        %
        [x,y,Omega,P,uMag,uX,uY,U,V,Fx,Fy,C] = import_Eulerian_Data(pathViz,numSim,Eulerian_Flags);


        % Imports Lagrangian Pt. FORCE (magnitude) DATA %
        %                      DEFINITIONS 
        %
        %      fX_Lag: forces in x-direction on boundary
        %      fY_Lag: forces in y-direction on boundary
        %       fLagMag: magnitude of force at boundary
        %   fLagNorm: magnitude of NORMAL force at boundary
        %   fLagTan: magnitude of TANGENT force at boundary
        %
        [fX_Lag,fY_Lag,fLagMag,fLagNorm,fLagTan] = import_Lagrangian_Force_Data(pathForce,numSim);
    
end