%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled  
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  nickabattista[at]gmail.com
% IB2d Created: May 27th, 2015
% Institution: TCNJ
%
% This code is capable of creating Lagrangian Structures using:
% 	1. Springs
% 	2. Beams (*torsional springs or non-invariant beams*)
% 	3. Target Points
%	4. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-Tension)")
%   .
%   .
%   .
%
% One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting lengths, etc
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
% 
% If you would like us to add a specific muscle model, please let Nick (nickabattista[at]gmail.co) know.
%
%--------------------------------------------------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: updates the beam attributes!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function beams_info = update_nonInv_Beams(dt,current_time,beams_info)


% beams_info:   col 1: 1ST PT.
%               col 2: MIDDLE PT. (where force is exerted)
%               col 3: 3RD PT.
%               col 4: beam stiffness
%               col 5: x-curvature
%               col 6: y-curvature

%IDs = beams_info(:,1);   % Gives Middle Pt.

% Time initialization
t1 = 0.25;                          % Period A->B (down-stroke)
t2 = 0.25;                          % Period B->A (up-stroke)
tR = 0.0;                           % Resting time between A->B and B->A
period = (t1+t2+tR);                % Time it takes to move to the right
t = rem(current_time,period);       % Recalculate Time (using modular arithmetic)

% Cubic Interpolation Information
p1 = 0.10;
p2 = 0.90;

% a COEFFICIENTS
a0 = 0; 
a1 = 0; 
a2 = 0; 
a3 = 11.111111111111088;

% b COEFFICIENTS
b0 = 0.013888888888889;
b1 = -0.416666666666666;
b2 = 4.166666666666658;
b3 = -2.777777777777770;

% c COEFFICIENTS
c0 =  -10.111111111111073;
c1 =  33.333333333333222;
c2 = -33.333333333333222;
c3 =  11.111111111111073;


% Read In y_Pts for two Phases!
[xP1,yP1,yP2] = read_File_In('swimmer.phases'); % NOTE xP1 = xP2
xP2 = xP1;

%
% FIRST WE COMPUTE THE INTERPOLATE GEOMETRY BETWEEN BOTH PHASES
%

%
% START THE INTERPOLATING BETWEEN STATES!
% 
if t <= t1 % STATE A -> STATE B

    % Scaling time for appropriate use in interp. function so tTilde\in[0,1]
    tTilde = (t/t1); 
    
    % Evaluate Pieceise Cubic Interpolation Poly
    if ( tTilde<=p1 )
        gFUNC = a0 + a1*tTilde + a2*tTilde^2 + a3*tTilde^3; 
    elseif ( (tTilde>p1) && (tTilde<=p2) )
        gFUNC = b0 + b1*tTilde + b2*tTilde^2 + b3*tTilde^3; 
    else
        gFUNC = c0 + c1*tTilde + c2*tTilde^2 + c3*tTilde^3; 
    end
    
    %xPts = xP1 + gFUNC*( xP2 - xP1 );	
    yPts = yP1 + gFUNC*( yP2 - yP1 );	
    
elseif ( t >= t1+tR ) % STATE B -> A
    
    % Scaling time for appropriate use in interp. function so tTilde\in[0,1]
    tTilde = (t-t1-tR)/(t2); 
    
    % Evaluate Pieceise Cubic Interpolation Poly
    if tTilde<=p1
        gFUNC = a0 + a1*tTilde + a2*tTilde^2 + a3*tTilde^3; 
    elseif tTilde<=p2
        gFUNC = b0 + b1*tTilde + b2*tTilde^2 + b3*tTilde^3; 
    else
        gFUNC = c0 + c1*tTilde + c2*tTilde^2 + c3*tTilde^3; 
    end
    
    %xPts = xP2 + gFUNC*( xP1 - xP2 );
    yPts = yP2 + gFUNC*( yP1 - yP2 );
    
end



%
% NOW WE UPDATE THE CURAVTURES APPROPRIATELY
%
%beams_info(:,5) = xPts(1:end-2)+xPts(3:end)-2*xPts(2:end-1);
beams_info(:,6) = yPts(1:end-2)+yPts(3:end)-2*yPts(2:end-1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: reads in info from file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x1,y1,y2] = read_File_In(file_name)

filename = file_name;  %Name of file to read in

fileID = fopen(filename);

    % Read in the file, use 'CollectOutput' to gather all similar data together
    % and 'CommentStyle' to to end and be able to skip lines in file.
    C = textscan(fileID,'%f %f %f','CollectOutput',1);

fclose(fileID);        %Close the data file.

mat_info = C{1};   %Stores all read in data

%Store all elements in matrix
mat = mat_info(1:end,1:end);

x1 =  mat(:,1); %store xVals1/2
y1 =  mat(:,2); %store yVals1 
y2 =  mat(:,3); %store yVals2