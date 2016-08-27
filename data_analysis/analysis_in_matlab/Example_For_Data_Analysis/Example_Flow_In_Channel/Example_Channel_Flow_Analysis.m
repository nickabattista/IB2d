%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
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
% FUNCTION: EXAMPLE DATA ANALYSIS CODE TO LOOK AT FLOW IN EMPTY CHANNEL
%
%     Note: 
%           (1) This code analyzes viz_IB2d code for channel flow in
%               /data_analysis/Example_For_Data_Analysis/Example_Flow_In_Channel/viz_IB2d
%           (2) Produces a plot of cross-sectional mag. of velocity for
%               different points along the channel at three times.
%           (3) USER-DEFINED functions are functions that users should make to
%               analyze their specific data sets
%           (4) MUST make sure to 'addpath' to where DA_Blackbox is, i.e.,
%               line 61
%           (5) MUST make sure to set path to desired dataset, i.e., in line
%               56
%          
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Example_Channel_Flow_Analysis()

% TEMPORAL INFO FROM input2d %
dt = 1e-4;      % Time-step
Tfinal = 0.015; % Final time in simulation
pDump=50;       % Note: 'print_dump' should match from input2d

% DATA ANALYSIS INFO %
start=1;                             % 1ST interval # included in data analysis
finish=3;                            % LAST interval # included in data analysis 
dump_Times = (start:1:finish)*pDump; % Time vector when data was printed in analysis

% SET PATH TO DESIRED viz_IB2d DATA %
path = 'viz_IB2d';

% SET PATH TO DA_BLACKBOX %
addpath('../../DA_Blackbox');

for i=start:1:finish
    
    % Points to desired data viz_IB2d data file
    if i<10
       numSim = ['000', num2str(i)];
    elseif i<100
       numSim = ['00', num2str(i)];
    elseif i<1000
       numSim = ['0', num2str(i)];
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
    %    Fx: x-directed Force     Fy: y-directed Force
    %
    %  Note: U(j,i): j-corresponds to y-index, i to the x-index
    %
    [x,y,Omega,P,uMag,uX,uY,U,V,Fx,Fy] = import_Eulerian_Data(path,numSim);

    %                                                                    %
    % *** USER DEFINED FUNCTIONS TO GET DESIRED ANALYSIS PT. INDICES *** %
    %                                                                    %
    if i==start
        xPts = [0.125 0.225 0.325 0.425];
        yPts = [0.405 0.595];
        [xInds,yInds] = give_Desired_Analysis_Points(x,y,xPts,yPts);
        vel_data = zeros(length(yInds),length(xInds),finish);
        %Inds = put_x_y_Indices_Together_For_Analysis(xInds,yInds);
    end
    
    %                                                                    %
    % ***** USER DEFINED FUNCTION TO SAVE DESIRED VELOCITY DATA *****    %
    %                                                                    %
    vel_data = store_Desired_Magnitude_Velocity_Data(uMag,vel_data,xInds,yInds,i);

end % END OF LOOP OVER ANALYZED TIME-PTS


%                                                                    %
% ***** USER DEFINED FUNCTION TO PLOT DESIRED VELOCITY DATA *****   %
%                                                                    %
yVals = y(yInds);
plot_Desired_Data(yVals,vel_data);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USER-FUNCTION: Finds desired analysis point indices
%
%        INPUTS: x: x-grid pts (row vector)
%                y: y-grid pts (column vector)
%                xPts: desired x-pts to analyze
%                yPts: desired y-pts to analyze
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xInds,yInds] = give_Desired_Analysis_Points(x,y,xPts,yPts)

% Get x-Indices
xInds = zeros(length(xPts),1);
for i=1:length(xPts)
   xPt = xPts(i);   % Get x-Pts
   k = 1;           % Initialize loop iteration variable
   while x(k) < xPt
      xInds(i) = k; 
      k = k+1;
   end
end

% Get y-Indices
yIndsAux = zeros(length(yPts),1);
for i=1:length(yPts)
   yPt = yPts(i);   % Get x-Pts
   k = 1;           % Initialize loop iteration variable
   while y(k) < yPt
      yIndsAux(i) = k; 
      k = k+1;
   end
end
yInds = yIndsAux(1):1:yIndsAux(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USER-FUNCTION: Rearranges desired analysis indices
%
%        INPUTS: xInds: indices for x-region of interest
%                yInds: indices for y-region of interest
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Inds = put_x_y_Indices_Together_For_Analysis(xInds,yInds)

% Order new index matrix %
Inds = zeros(length(xInds)*length(yInds),2);
count = 1;
for i=1:length(xInds)
    for j=1:length(yInds)
        Inds(count,1) = xInds(i);
        Inds(count,2) = yInds(j);
        count=count+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USER-FUNCTION: Stores desired magnitude of velocity data
%
%        INPUTS: uMag: magnitude of velocity from simulation
%                vel_data: stored mag. of velocity data in 3D-matrix
%                xInds/yInds: indices of where to save data
%                i: ith time-step to store data from
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vel_data = store_Desired_Magnitude_Velocity_Data(uMag,vel_data,xInds,yInds,i)

% NOTE: vel_data(k,i): : k: magnitude of velocity data  (row)
%                      : j: xPt you're storing data for (column)
%                      : i: ith level you're storing in (level)

for j=1:length(xInds)
    for k=1:length(yInds)
       vel_data(k,j,i) = uMag( yInds(k),xInds(j) ); 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USER-FUNCTION: Plots magnitude of velocity data
%
%        INPUTS: vel_data: stored mag. of velocity data in 3D-matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_Desired_Data(yVals,vel_data)

fs = 16; %Font Size

% Set Figure Size
FigHand = figure(1);
set(FigHand,'Position',[100,100,1024,895]);
%
% Make Figure!
%
subplot(3,1,1)
mat = vel_data(:,:,1);
maxVal = max(max(mat));
plot(yVals,vel_data(:,1,1),'*-'); hold on;
plot(yVals,vel_data(:,2,1),'r*-'); hold on;
plot(yVals,vel_data(:,3,1),'g*-'); hold on;
plot(yVals,vel_data(:,4,1),'k*-'); hold on;
axis([0.4 0.6 0 1.1*maxVal]);
leg=legend('x=0.125','x=0.175','x=0.225','x=0.275');
title('t=0.005'); 
ylabel('Mag. Velocity','FontSize',fs); xlabel('y','FontSize',fs);
set(leg,'FontSize',fs);
set(gca,'FontSize',fs-1);
%
subplot(3,1,2)
mat = vel_data(:,:,2);
maxVal = max(max(mat));
plot(yVals,vel_data(:,1,2),'*-');  hold on;
plot(yVals,vel_data(:,2,2),'r*-'); hold on;
plot(yVals,vel_data(:,3,2),'g*-'); hold on;
plot(yVals,vel_data(:,4,2),'k*-'); hold on;
axis([0.4 0.6 0 1.1*maxVal]);
leg=legend('x=0.125','x=0.175','x=0.225','x=0.275');
ylabel('Mag. Velocity','FontSize',fs); xlabel('y','FontSize',fs);
set(leg,'FontSize',fs);
set(gca,'FontSize',fs-1);
title('t=0.01');
%
subplot(3,1,3)
mat = vel_data(:,:,3);
maxVal = max(max(mat));
plot(yVals,vel_data(:,1,3),'*-'); hold on;
plot(yVals,vel_data(:,2,3),'r*-'); hold on;
plot(yVals,vel_data(:,3,3),'g*-'); hold on;
plot(yVals,vel_data(:,4,3),'k*-'); hold on;
axis([0.4 0.6 0 1.1*maxVal]);
leg=legend('x=0.125','x=0.175','x=0.225','x=0.275');
ylabel('Mag. Velocity','FontSize',fs); xlabel('y','FontSize',fs);
set(leg,'FontSize',fs);
set(gca,'FontSize',fs-1);
title('t=0.015');



