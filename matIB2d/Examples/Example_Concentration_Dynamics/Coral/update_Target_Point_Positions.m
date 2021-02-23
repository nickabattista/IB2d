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
%	4. Muscle-Model (combined Force-Length-Velocity model, "Hill+(Length-Tension)")
%
% One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting lengths, etc
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
% 
% If you would like us to add a specific muscle model, please let Nick (nick.battista@unc.edu) know.
%
%--------------------------------------------------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: updates the target point positions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function targets = update_Target_Point_Positions(dt,current_time,targets)

dx=2/256;

xoffset = 1; % To put geometry into QUADRANT-1
yoffset = 0.5; % To put geometry into QUADRANT-1

% Note: (1) ds for stored data is 0.6/(2*1024)
%       (2) 'ratio' is comparing 1024:desired resolution
ds=dx/2;


load('total_coeffs.mat')
load('coral_coeff_30.mat')

%arclength
s=(0:ds:total_meanL)/total_meanL;

%IDs = targets(:,1);                 % Stores Lag-Pt IDs in col vector
%xPts= targets(:,2);                % Original x-Values of x-Target Pts.
%yPts= targets(:,3);                 % Original y-Values of y-Target Pts.
%kStiffs = targets(:,4);             % Stores Target Stiffnesses 

indx=ceil(current_time/dt)+1;

load('cval.mat')

C1=c1vals(indx,:);
C2=c2vals(indx,:);

     XbL_1=(C1(1)*s.^3+C1(2)*s.^2+C1(3)*s+C1(4));
     YbL_1=(C2(1)*s.^3+C2(2)*s.^2+C2(3)*s+C2(4));

     L_ten=sum(sqrt((XbL_1(2:end)-XbL_1(1:end-1)).^2 +(YbL_1(2:end)-YbL_1(1:end-1)).^2 ));

     XbL=(C1(1)*s.^3+C1(2)*s.^2+C1(3)*s+C1(4))*total_meanL/L_ten+total_offset;
     YbL=(C2(1)*s.^3+C2(2)*s.^2+C2(3)*s+C2(4))*total_meanL/L_ten;

    XbR=-XbL;
    YbR=YbL;


% XbStem=(XbR(1)+ds):ds:(XbL(1)-ds);
YbStem=(YbL(1))+0*XbStem;

%length(XbStem)
% XbStem=linspace(XbR(1)+ds,XbL(1)-ds,Nstem);
% YbStem=(YbL(1))+0*XbStem;


x1=[flip(XbR) XbStem XbL];
y1=[flip(YbR) YbStem YbL];

% Put Geometry Together for One Polyp
x = x1+xoffset;
y = y1+yoffset;



    
    targets(:,2) = x;      % Store new xVals
    targets(:,3) = y; % Store new yVals
end

