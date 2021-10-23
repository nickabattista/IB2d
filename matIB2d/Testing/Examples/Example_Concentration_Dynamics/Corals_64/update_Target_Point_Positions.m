%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  battistn@tcnj.edu
% Date Modified: April 2021
% Current Institution: TCNJ
%
% IB2d Date Created: May 27th, 2015
% Institution Created: UNC-CH
%
% This code is capable of creating Lagrangian Structures using:
% 	1. Springs
% 	2. Beams (torsional springs or non-invariant beams)
% 	3. Target Points
%	4. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-Tension)")
%
% One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting lengths, etc
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
% 
% If you would like us to add a specific muscle model, please let Nick (battistn@tcnj.edu) know.
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: updates the target point positions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function targets = update_Target_Point_Positions(dt,current_time,targets)

%IDs = targets(:,1);                 % Stores Lag-Pt IDs in col vector
%xPts= targets(:,2);                % Original x-Values of x-Target Pts.
%yPts= targets(:,3);                 % Original y-Values of y-Target Pts.
%kStiffs = targets(:,4);             % Stores Target Stiffnesses 

%-----------------------------------------------------
% Geometric Translations to Quadrant 1
%    since stored .mat geo data centered at (0,0)
%-----------------------------------------------------
xoffset = 1;   % To put geometry into QUADRANT-1
yoffset = 0.5; % To put geometry into QUADRANT-1

%-----------------------------------------------
% Note: (1) ds for stored data is 0.6/(2*1024)
%       (2) 'ratio' is comparing 1024:desired resolution
%-----------------------------------------------
Lx=2;          % Horizontal Eulerian grid length
Nx=64;        % Eulerian grid resolution
dx=Lx/Nx;      % Eulerian grid spacing
ds=dx/2;       % Lagrangian point spacing


%-----------------------------------------------
% Load prescribed position data for tentacles
%-----------------------------------------------
load('total_coeffs.mat')
load('coral_coeff_30.mat')

%-----------------------------------------------
% Arclength
%-----------------------------------------------
s=(0:ds:total_meanL)/total_meanL;

%---------------------------------------------------------
% Get index correponding to current time in 
%   simulation to determine how interpolation occurs
%---------------------------------------------------------
indx=ceil(current_time/dt)+1;

%-----------------------------------------------
% Load geometry state data
%-----------------------------------------------
load('cval.mat')

%-------------------------------------------------
% Get interpolation polynomial coefficients and
%   then determine geometry of LEFT tentacle
%-------------------------------------------------
C1=c1vals(indx,:);
C2=c2vals(indx,:);

XbL_1=(C1(1)*s.^3+C1(2)*s.^2+C1(3)*s+C1(4));
YbL_1=(C2(1)*s.^3+C2(2)*s.^2+C2(3)*s+C2(4));
L_ten=sum(sqrt((XbL_1(2:end)-XbL_1(1:end-1)).^2 +(YbL_1(2:end)-YbL_1(1:end-1)).^2 ));

XbL=(C1(1)*s.^3+C1(2)*s.^2+C1(3)*s+C1(4))*total_meanL/L_ten+total_offset;
YbL=(C2(1)*s.^3+C2(2)*s.^2+C2(3)*s+C2(4))*total_meanL/L_ten;

%-------------------------------------------------------
% Set up symmetric RIGHT tentacle interpolation states
%-------------------------------------------------------
XbR=-XbL;
YbR=YbL;

%-------------------------------------------------------
% Get coral polyp STEM geometry
%-------------------------------------------------------
YbStem=(YbL(1))+0*XbStem;

%-------------------------------------------------------
% Combine geometry into one vector and translate
%-------------------------------------------------------
x=[flip(XbR) XbStem XbL];
y=[flip(YbR) YbStem YbL];
%
x = x+xoffset;
y = y+yoffset;

%-------------------------------------------------------
% Update target point positions!
%-------------------------------------------------------
targets(:,2) = x;   % Store new xVals
targets(:,3) = y;   % Store new yVals


