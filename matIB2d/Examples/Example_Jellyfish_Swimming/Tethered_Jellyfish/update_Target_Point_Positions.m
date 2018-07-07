%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
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
%	4. Muscle-Model (combined Force-Length-Velocity model, "Hill+(Length-Tension)")
%
% One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting lengths, etc
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
% 
% If you would like us to add a specific muscle model, please let Nick (nickabattista@gmail.com) know.
%
%--------------------------------------------------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: updates the target point positions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function targets = update_Target_Point_Positions(dt,current_time,targets)


%IDs = targets(:,1);                 % Stores Lag-Pt IDs in col vector
%xPts= targets(:,2);                 % Original x-Values of x-Target Pts.
%yPts= targets(:,3);                 % Original y-Values of y-Target Pts.
%kStiffs = targets(:,4);             % Stores Target Stiffnesses 
%N_target = length(targets(:,1));    % Gives total number of target pts!


tR1 = 0.25;                          % Period A->B (and B->A) for RIGHT
tR2 = 0.25;                          % Period B->C (abd C->B) for RIGHT
periodR = (tR1+tR2);                % Time it takes to move to the right

off_Left = 0;%0.10*periodR;             % OFFSET FOR LEFT ARM
tL1 = 0.25;                          % Period A->B (and B->A) for LEFT
tL2 = 0.25;                          % Period B->C (abd C->B) for LEFT
periodL = (tL1+tL2);        % Time it takes to move to the LEFT

tR = rem(current_time,periodR);      % "Time" (moded out) according to right arm

tL = rem(current_time,periodL);      % "Time" (moded out) according to left arm


% Cubic Interpolation Information
p1 = 0.125;
p2 = 0.875;

% a COEFFICIENTS
a0 = 0; 
a1 = 0; 
a2 = 0; 
a3 = 9.1429;

% b COEFFICIENTS
b0 = 0.0238;
b1 = -0.5714;
b2 = 4.5714;
b3 = -3.0476;

% c COEFFICIENTS
c0 = -8.1429;
c1 = 27.4286;
c2 = -27.4286;
c3 = 9.1429;

%
% JELLYFISH STATES: 
%   XY(:,1) = 1st X-positions
%   XY(2:150,2) = 2nd X-positions (RIGHT)
%   XY(151:end,2) = 2nd X-positions (LEFT)
%
%   XY(:,3) = 1st Y-positions
%   XY(2:150,4) = 2nd Y-positions (RIGHT)
%   XY(151:end,4) = 2nd Y-positions (LEFT)
%
XY = read_In_Phase_Positions('XY_2Pos.txt');

% Coming from give_Me_Jelly:
NArm = 174;

% 1st and 2nd position of RIGHT arm
xR1 = XY(2:NArm+1,1); yR1 = XY(2:NArm+1,3);
xR2 = XY(2:NArm+1,2); yR2 = XY(2:NArm+1,4);

% 1st and 2nd position of LEFT arm
xL1 = XY(NArm+2:end,1); yL1 = XY(NArm+2:end,3);
xL2 = XY(NArm+2:end,2); yL2 = XY(NArm+2:end,4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START THE INTERPOLATING BETWEEN STATES!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%
%
% RIGHT ARM!
%
%%%%%%%%%%%%%

if tR <= tR1+off_Left % STATE A -> STATE B
    
    % Scaling time for appropriate use in interp. function so tRTilde\in[0,1]
    tRTilde = ( tR / tR1 ); 
    
    % Evaluate Piecewise Cubic Interpolation Poly
    if tRTilde<=p1
        gFUNC = a0 + a1*tRTilde + a2*tRTilde^2 + a3*tRTilde^3; 
    elseif tRTilde<=p2
        gFUNC = b0 + b1*tRTilde + b2*tRTilde^2 + b3*tRTilde^3; 
    else
        gFUNC = c0 + c1*tRTilde + c2*tRTilde^2 + c3*tRTilde^3; 
    end
    
    % MOVING X-POSITIONS
    targets(2:NArm+1,2) = xR1 + gFUNC*( xR2 - xR1 );
    
    % MOVING Y-POSITIONS
    targets(2:NArm+1,3) = yR1 + gFUNC*( yR2 - yR1 );
        
elseif tR <= (tR1+tR2) % STATE B -> A
    
    % Scaling time for appropriate use in interp. function so tRTilde\in[0,1]
    tRTilde = (tR-tR1)/( tR2 ); 
    
    % Evaluate Piecewise Cubic Interpolation Poly
    if tRTilde<=p1
        gFUNC = a0 + a1*tRTilde + a2*tRTilde^2 + a3*tRTilde^3; 
    elseif tRTilde<=p2
        gFUNC = b0 + b1*tRTilde + b2*tRTilde^2 + b3*tRTilde^3; 
    else
        gFUNC = c0 + c1*tRTilde + c2*tRTilde^2 + c3*tRTilde^3; 
    end

    % MOVING X-POSITIONS
    targets(2:NArm+1,2) = xR2 + gFUNC*( xR1 - xR2 );
    
    % MOVING Y-POSITIONS
    targets(2:NArm+1,3) = yR2 + gFUNC*( yR1 - yR2 );
    
end


%%%%%%%%%%%%%
%
% LEFT ARM!
%
%%%%%%%%%%%%%

if current_time >= off_Left

    if tL <= tL1 % STATE A -> STATE B

        % Scaling time for appropriate use in interp. function so tRTilde\in[0,1]
        tLTilde = (tL/tL1); 

        % Evaluate Piecewise Cubic Interpolation Poly
        if tLTilde<=p1
            gFUNC = a0 + a1*tLTilde + a2*tLTilde^2 + a3*tLTilde^3; 
        elseif tRTilde<=p2
            gFUNC = b0 + b1*tLTilde + b2*tLTilde^2 + b3*tLTilde^3; 
        else
            gFUNC = c0 + c1*tLTilde + c2*tLTilde^2 + c3*tLTilde^3; 
        end

        % MOVING X-POSITIONS
        targets(NArm+2:end,2) = xL1 + gFUNC*( xL2 - xL1 );

        % MOVING Y-POSITIONS
        targets(NArm+2:end,3) = yL1 + gFUNC*( yL2 - yL1 );

    elseif tL <= (tL1+tL2) % STATE B -> A

        % Scaling time for appropriate use in interp. function so tRTilde\in[0,1]
        tLTilde = (tL-tL1)/(tL2); 

        % Evaluate Piecewise Cubic Interpolation Poly
        if tLTilde<=p1
            gFUNC = a0 + a1*tLTilde + a2*tLTilde^2 + a3*tLTilde^3; 
        elseif tLTilde<=p2
            gFUNC = b0 + b1*tLTilde + b2*tLTilde^2 + b3*tLTilde^3; 
        else
            gFUNC = c0 + c1*tLTilde + c2*tLTilde^2 + c3*tLTilde^3; 
        end

        % MOVING X-POSITIONS
        targets(NArm+2:end,2) = xL2 + gFUNC*( xL1 - xL2 );

        % MOVING Y-POSITIONS
        targets(NArm+2:end,3) = yL2 + gFUNC*( yL1 - yL2 );

    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in the # of vertex pts and all the vertex pts from the
%           .vertex file.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PTS = read_In_Phase_Positions(struct_name)


filename = struct_name;  %Name of file to read in
fileID = fopen(filename);

% Read in the file, use 'CollectOutput' to gather all similar data together
% and 'CommentStyle' to to end and be able to skip lines in file.
C = textscan(fileID,'%f %f %f %f','CollectOutput',1);

fclose(fileID);     %Close the data file.

vertices = C{1};    %Stores all read in data in vertices (N+1,2) array

PTS = vertices(1:end,1:4);
