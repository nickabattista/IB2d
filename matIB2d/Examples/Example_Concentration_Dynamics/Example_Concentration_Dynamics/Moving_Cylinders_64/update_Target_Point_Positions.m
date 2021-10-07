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
%xPts= targets(:,2);                 % Original x-Values of x-Target Pts.
%yPts= targets(:,3);                 % Original y-Values of y-Target Pts.
%kStiffs = targets(:,4);             % Stores Target Stiffnesses 
%N_target = length(targets(:,1));    % Gives total number of target pts!


% Time Information
t1 = 0.5;                            % Period A -> B && B -> A                          
double_period = 2*t1;                % Time it takes for A -> B -> A
t = rem(current_time,double_period); % Time (moded out)


% Cubic Interpolation Information
p1 = 0.05;
p2 = 0.95;

% a COEFFICIENTS
a0 = 0;
a1 = 0;
a2 = 0;
a3 = 21.052631578947341;

% b COEFFICIENTS
b0 = 0.002923976608187;
b1 = -0.175438596491228;
b2 = 3.508771929824557;
b3 = -2.339181286549700;

% c COEFFICIENTS
c0 = -20.052631578947224;
c1 =  63.157894736841683;
c2 = -63.157894736841683;
c3 = 21.052631578947228;

A = read_In_State('State_A.pts');
B = read_In_State('State_B.pts');

%
% START THE INTERPOLATING BETWEEN STATES!
% 
if t <= t1 % STATE A -> STATE B
    
    % Scaling time for appropriate use in interp. function so tTilde\in[0,1]
    tTilde = (t/t1); 
    
    % Evaluate Pieceise Cubic Interpolation Poly
    if tTilde<=p1
        gFUNC = a0 + a1*tTilde + a2*tTilde^2 + a3*tTilde^3; 
    elseif tTilde<=p1+p2
        gFUNC = b0 + b1*tTilde + b2*tTilde^2 + b3*tTilde^3; 
    else
        gFUNC = c0 + c1*tTilde + c2*tTilde^2 + c3*tTilde^3; 
    end
    
    targets(:,2) = A(:,1) + gFUNC*( B(:,1) - A(:,1) );
    targets(:,3) = A(:,2) + gFUNC*( B(:,2) - A(:,2) );
    
else % STATE B -> A
   
    A = read_In_State('State_A.pts');
    B = read_In_State('State_B.pts');
    
    % Scaling time for appropriate use in interp. function so tTilde\in[0,1]
    tTilde = (t-t1)/(t1); 
    
    % Evaluate Pieceise Cubic Interpolation Poly
    if tTilde<=p1
        gFUNC = a0 + a1*tTilde + a2*tTilde^2 + a3*tTilde^3; 
    elseif tTilde<=p1+p2
        gFUNC = b0 + b1*tTilde + b2*tTilde^2 + b3*tTilde^3; 
    else
        gFUNC = c0 + c1*tTilde + c2*tTilde^2 + c3*tTilde^3; 
    end
    
    targets(:,2) = B(:,1) + gFUNC * ( A(:,1) - B(:,1) );
    targets(:,3) = B(:,2) + gFUNC * ( A(:,2) - B(:,2) );
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in the # of vertex pts and all the vertex pts from the
%           .vertex file.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PTS = read_In_State(struct_name)


filename = struct_name;  %Name of file to read in
fileID = fopen(filename);

% Read in the file, use 'CollectOutput' to gather all similar data together
% and 'CommentStyle' to to end and be able to skip lines in file.
C = textscan(fileID,'%f %f','CollectOutput',1);

fclose(fileID);     %Close the data file.

vertices = C{1};    %Stores all read in data in vertices (N+1,2) array

PTS = vertices(1:end,1:2);
