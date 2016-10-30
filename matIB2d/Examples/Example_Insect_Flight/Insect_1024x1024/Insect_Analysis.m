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


function Insect_Analysis()

% TEMPORAL INFO FROM input2d %
dt = 2.5e-5;      % Time-step
Tfinal = 0.05;   % Final time in simulation
pDump=40;        % Note: 'print_dump' should match from input2d

% DATA ANALYSIS INFO %
start=1;                             % 1ST interval # included in data analysis
finish=25;                           % LAST interval # included in data analysis 
dump_Times = (start:1:finish)*pDump; % Time vector when data was printed in analysis

% SET PATH TO DESIRED viz_IB2d DATA %
pathMain = '/Users/nick_battista/Desktop/IB2d/matIB2d/Examples/Example_Insect_Flight/Insect_1024x1024';

% SET PATH TO DA_BLACKBOX %
addpath('/Users/nick_battista/Desktop/IB2d/data_analysis/analysis_in_matlab/DA_Blackbox');

    
pathViz = [pathMain '/viz_IB2d/'];
pathForce = [pathMain '/hier_IB2d_data'];

%simVec = {'32x32','64x64','96x96','128x128','256x256','512x512','768x768','1024x1024'};

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
    %[xLag,yLag] = give_Lag_Positions(pathViz,numSim);

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
    %[x,y,Omega,P,uMag,uX,uY,U,V,Fx,Fy] = import_Eulerian_Data(pathViz,numSim);


    % Imports Lagrangian Pt. FORCE (magnitude) DATA %
    %                      DEFINITIONS 
    %
    %      fX_Lag: forces in x-direction on boundary
    %      fY_Lag: forces in y-direction on boundary
    %       fLagMag: magnitude of force at boundary
    %   fLagNorm: magnitude of NORMAL force at boundary
    %   fLagTan: magnitude of TANGENT force at boundary
    %
    [fX_Lag,fY_Lag] = import_Lagrangian_Force_Data_Insect(pathForce,numSim);


    fx_Mat(i,1) = mean( fX_Lag );
    fy_Mat(i,1) = mean( fY_Lag );

end

cd ..;

strName = 'fX_1024';
print_Matrix_To_Txt_File(fx_Mat,strName)

strName = 'fY_1024';
print_Matrix_To_Txt_File(fy_Mat,strName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION prints matrix to file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Matrix_To_Txt_File(a,strName)

nameTxt = [strName '.txt'];

fid = fopen(nameTxt, 'wt'); % Open for writing
for i=1:size(a,1)
   fprintf(fid, '%d ', a(i,:));
   fprintf(fid, '\n');
end
fclose(fid);
    