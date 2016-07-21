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
%               line 56
%           (3) MUST make sure to set path to desired dataset, i.e., in line
%               59
%          
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Standard_Data_Analysis_Script()

% TEMPORAL INFO FROM input2d %
dt = 5e-5;      % Time-step
Tfinal = 0.1;   % Final time in simulation
pDump=50;       % Note: 'print_dump' should match from input2d

% DATA ANALYSIS INFO %
start=1;                             % 1ST interval # included in data analysis
finish=40;                           % LAST interval # included in data analysis 
dump_Times = (start:1:finish)*pDump; % Time vector when data was printed in analysis

% SET PATH TO DESIRED viz_IB2d DATA %
path = '/Users/nick_battista/Desktop/IB2d/matIB2d/Examples/Example_Porous_Rubberbands/Compare_Porous_NonPorous_Rubberbands/viz_IB2d';

% SET PATH TO DA_BLACKBOX %
addpath('../DA_Blackbox');

for i=start:1:finish
    
    % Points to desired data viz_IB2d data file
    if i<10
       numSim = ['000', num2str(i)];
    elseif i<100
       numSim = ['00', num2str(i)];
    elseif i<1000
       numSim = ['0', num2str(i) '.vtk'];
    else
       numSim = num2str(i);
    end
    
    % Imports immersed boundary positions %
    [xLag,yLag] = give_Lag_Positions(path,numSim);

    % Imports (x,y) grid values and ALL Eulerian Data %
    %                      DEFINITIONS 
    %          x: x-grid                y: y-grid
    %       Omega: vorticity           P: pressure
    %    uMag: mag. of velocity  
    %    uX: mag. of x-Velocity   uY: mag. of y-Velocity  
    %    U: x-directed velocity   V: y-directed velocity
    % 
    %  Note: U(j,i): j-corresponds to y-index, i to the x-index
    %
    [x,y,Omega,P,uMag,uX,uY,U,V] = import_Eulerian_Data(path,numSim);

end