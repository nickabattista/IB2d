%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled  
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  nickabattista@gmail.com
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
% If you would like us to add a specific muscle model, please let Nick (nickabattista@gmail.com) know.
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

function Sawtooth_Data_Analysis_AvgEverything()

% TEMPORAL INFO FROM input2d %
dt = 2.5e-5;     % Time-step
Tfinal = 10.0;   % Final time in simulation
pDump=400;       % Note: 'print_dump' should match from input2d

% DATA ANALYSIS INFO %
start=1;                             % 1ST interval # included in data analysis
finish=999;                         % LAST interval # included in data analysis 
dump_Times = (start:1:finish)*pDump; % Time vector when data was printed in analysis

% INITIALIZE STORAGE MATRICES
flow_MAT = zeros(finish,1);
vol_MAT = flow_MAT;

% VARIABLES TO NON-DIMENSIONALIZE
G = 0.325; % cm

% SET PATH TO DESIRED viz_IB2d / hier_IB2d_data DATA %
pathViz = 'viz_IB2d';
%pathForce= 'hier_IB2d_data';

% SET PATH TO DA_BLACKBOX %
addpath('/Users/battistn/Desktop/IB2d/data_analysis/analysis_in_matlab/DA_Blackbox');  % <--- MAKE SURE THIS IS PASSED ACCORDINGLY

%
% START LOOP FOR DATA ANALYSIS OF SPECIFIC SIMULATION
%
fprintf('\n\nStartin Data Analysis...\n\n');
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

        if mod(i,10) == 0
            fprintf('Analyzing %d of %d\n',i,finish);
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
        Eulerian_Flags(1) = 1;   % OMEGA
        Eulerian_Flags(2) = 0;   % PRESSURE
        Eulerian_Flags(3) = 0;   % uMAG
        Eulerian_Flags(4) = 0;   % uX (mag. x-component of velocity)
        Eulerian_Flags(5) = 0;   % uY (mag. x-component of velocity)
        Eulerian_Flags(6) = 1;   % uVEC (vector components of velocity: U,V)
        Eulerian_Flags(7) = 0;   % Fx (x-component of force )
        Eulerian_Flags(8) = 0;   % Fy (y-component of force)
        %
        %[x,y,Omega,P,uMag,uX,uY,U,V,Fx,Fy] = import_Eulerian_Data(pathViz,numSim,Eulerian_Flags);
        [x,y,~,~,~,~,~,U,~,~,~] = import_Eulerian_Data(pathViz,numSim,Eulerian_Flags);


        % Imports Lagrangian Pt. FORCE (magnitude) DATA %
        %                      DEFINITIONS 
        %
        %      fX_Lag: forces in x-direction on boundary
        %      fY_Lag: forces in y-direction on boundary
        %       fLagMag: magnitude of force at boundary
        %   fLagNorm: magnitude of NORMAL force at boundary
        %   fLagTan: magnitude of TANGENT force at boundary
        %
        %[fX_Lag,fY_Lag,fLagMag,fLagNorm,fLagTan] = import_Lagrangian_Force_Data(pathForce,numSim);
    
        % give y Value of top-rightmost point of sawtooth
        yTop = yLag(end/2);
        yBot = yLag(end);
        xRight = xLag(end/2);
        xLeft = xLag(1);
        
        % gives inds = [xLeft xRight yBot yTop] <-- all indices
        inds = give_Correct_Inds(x,y,xLeft,xRight,yBot,yTop);
        
        % define any buffers, if desired
        xBuff = 3; % how wide of region to average over ("xBuff-boxes wide")
        yBuff = 0; % how vertical of region to average over
        
        % gives vertical distance to average over
        vert_dist = y(inds(4)-yBuff)-y(inds(3)+yBuff);
        
        % get inds for averaging
        inds_Vert = give_Inds_For_U_V(inds,xBuff,yBuff);
        
        % get flows left/right
        flow = get_Average_Flows_Left_Right(inds_Vert,inds,U);
        
        % flow (speed) matrix
        flow_MAT(i) = flow;
        %flow_MAT(i,2) = flow_R;
        
        % volumetric flow matrix ( velocity x width of channel at that time)
        vol_MAT(i) = flow * vert_dist / G;
        %vol_MAT(i,2) = flow_R * vert_dist / G;
        
end



flow_MAT = [ 0; flow_MAT];

vol_MAT = [0; vol_MAT];

mean( flow_MAT(2:end,1) )
%mean( flow_MAT(2:end,2) )

mean( vol_MAT(2:end,1) )
%mean( vol_MAT(2:end,2) )

% PRINT Information for Simulation %
% strName = ['FLOW_DATA.txt'];
% fileID = fopen(strName,'w');
% for i=1:finish
%     fprintf(fileID,'%1.16e %1.16e %1.16 %1.16e\n',flow_MAT(i,1),flow_MAT(i,2),vol_MAT(i,1),vol_MAT(i,2) );
% end
% fclose(fileID);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: averages flows over region
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function flow = get_Average_Flows_Left_Right(inds_vert,inds,U)

len = length(inds_vert);

mean_tot = 0;
for i=1:len
    for j=inds(3):inds(4)
        mean_tot = mean_tot + U( inds_vert(i), j );
    end
end

flow = mean_tot / ( len*(inds(4)-inds(3)+1) );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives (xInd,yInd) for data points
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function inds_Vert = give_Inds_For_U_V(inds,xBuff,yBuff)

    inds_Vert = [inds(3)+yBuff:1:inds(4)-yBuff]';
%     lenV = length(inds_Vert);
%     inds_hor = inds(1)*ones(lenV,1);
%     inds_horR= inds(2)*ones(lenV,1);
%     
%     for i=2:xBuff
%         inds_hor = [inds_hor; (inds(1)+(i-1))*ones(lenV,1)];
%         inds_horR = [inds_horR; (inds(2)-(i-1))*ones(lenV,1)];
%     end
%     
%     inds_Vert = [inds_Vert; inds_Vert; inds_Vert];
% 
%     % Combine: yInds are first column, xInds are second and third column
%     inds_UV = [inds_Vert inds_hor inds_horR];
% 
%     clear inds_hor inds_horR inds_Vert

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives correct inds to sweep for averaging
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function inds = give_Correct_Inds(x,y,xLeft,xRight,yBot,yTop)


%
% Find correct lower y-Ind
%
y1 = y(1);
j=1;
while y1 < yBot
   j=j+1;
   y1 = y(j);
end

indY_bot = j;

%
% Find correct upper y-Ind
%
y1 = y(end);
j=length(y);
while y1 > yTop
   j=j-1;
   y1 = y(j);
end

indY_top = j;


%
% Find correct left x-Ind
%
x1 = x(1);
j=1;
while x1 < xLeft
   j=j+1;
   x1 = x(j);
end

indX_left = j;


%
% Find correct RIGHT x-Ind
%
x1 = x(end);
j=length(x);
while x1 > xRight
   j=j-1;
   x1 = x(j);
end

indX_right = j;


inds = [indX_left indX_right indY_bot indY_top];