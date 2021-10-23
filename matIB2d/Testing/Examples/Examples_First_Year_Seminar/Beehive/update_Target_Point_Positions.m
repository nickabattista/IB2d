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

% INFO FROM TARGET POINTS %
IDs = targets(:,1);                 % Stores Lag-Pt IDs in col vector
xPts= targets(:,2);                 % Original x-Values of x-Target Pts.
yPts= targets(:,3);                 % Original y-Values of y-Target Pts.
kStiffs = targets(:,4);             % Stores Target Stiffnesses 


% -> FINDING NEXT ANGULAR FREQUENCY <- % 
[angInput1,angInput2] = give_Me_Oscillatory_Angle(dt,current_time);

% Declare how open you want the movements to go
full_ang = pi/3;
dbl_ang = 2*full_ang;

% Find associated angle
ang1 = mod( angInput1, dbl_ang );
if ( (ang1 > full_ang ) && (ang1 <= dbl_ang) )
    ang1 = dbl_ang-ang1;
end

ang2 = mod( angInput2, dbl_ang );
if ( (ang2 > full_ang ) && (ang2 <= dbl_ang) )
    ang2 = dbl_ang-ang2;
end


% Read In REFERENCE Pts!
[xRef,yRef] = read_File_In('hive.vertex');

% Store Values for Centers of Rotation
xL_1 = xRef(1);       yL_1 = yRef(1);           % Left side of Left Pair
xR_1 = xRef(end/4+1); yR_1 = yRef(end/4+1);     % Right side of Left Pair

xL_2 = xRef(end/2+1);   yL_2 = yRef(end/2+1);   % Left side of Right Pair
xR_2 = xRef(3*end/4+1); yR_2 = yRef(3*end/4+1); % Right side of Right Pair

% -> Rotate Geometry <- %
% LEFT PAIR %
[xR_Ref_1,yR_Ref_1] = rotate_Geometry(-ang1,xR_1,yR_1,xRef(end/4+1:end/2),yRef(end/4+1:end/2) );
[xL_Ref_1,yL_Ref_1] = rotate_Geometry(ang1,xL_1,yL_1,xRef(1:end/4),yRef(1:end/4) );
% RIGHT PAIR %
[xR_Ref_2,yR_Ref_2] = rotate_Geometry(-ang2,xR_2,yR_2,xRef(3*end/4+1:end),yRef(3*end/4+1:end) );
[xL_Ref_2,yL_Ref_2] = rotate_Geometry(ang2,xL_2,yL_2,xRef(end/2+1:3*end/4),yRef(end/2+1:3*end/4) );

% Store New Geometry
targets(IDs,2) = [xL_Ref_1 xR_Ref_1 xL_Ref_2 xR_Ref_2]; 
targets(IDs,3) = [yL_Ref_1 yR_Ref_1 yL_Ref_2 yR_Ref_2]; 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: rotate geometry
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,y] = rotate_Geometry(ang,xC,yC,xRef,yRef)

len = length(xRef);
x = zeros(1,len); y=x;

xRef = xRef - xC;
yRef = yRef - yC;

for i=1:len
   x(i) =  xRef(i)*cos(ang) - yRef(i)*sin(ang);
   y(i) =  xRef(i)*sin(ang) + yRef(i)*cos(ang);
end

x = x + xC;
y = y + yC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: reads in info from file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x1,y1] = read_File_In(file_name)

filename = file_name;  %Name of file to read in

fileID = fopen(filename);

    % Read in the file, use 'CollectOutput' to gather all similar data together
    % and 'CommentStyle' to to end and be able to skip lines in file.
    C = textscan(fileID,'%f %f','CollectOutput',1);

fclose(fileID);        %Close the data file.

mat_info = C{1};   %Stores all read in data

%Store all elements in matrix
mat = mat_info(2:end,1:end);

x1 =  mat(:,1); %store xVals 
y1 =  mat(:,2); %store yVals


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives me next angle
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [angInput1,angInput2] = give_Me_Oscillatory_Angle(dt,current_time)

n = ceil((current_time+dt) / dt);

dtheta = pi*dt;
angVec1 = dtheta:dtheta:2*pi;
n1 = mod(n,length(angVec1));

dtheta = pi*dt/2;
angVec2 = dtheta:dtheta:2*pi;
n2 = mod(n,length(angVec2));


angInput1 = angVec1(n1);
angInput2 = angVec2(n2);