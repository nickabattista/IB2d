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
% FUNCTION: updates the target point positions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function targets = update_Target_Point_Positions(dt,current_time,targets)


%IDs = targets(:,1);                 % Stores Lag-Pt IDs in col vector
%xPts= targets(:,2);                 % Original x-Values of x-Target Pts.
%yPts= targets(:,3);                 % Original y-Values of y-Target Pts.
%kStiffs = targets(:,4);             % Stores Target Stiffnesses 
%N_target = length(targets(:,1));    % Gives total number of target pts!


t1 = 0.01;
t2 = 0.02;
period = 2*(t1+t2);                 % Time it takes to move to the right
t = rem(current_time,period);       % Time (moded out)


if t <= t1
    A = read_In_State('State_A.pts');
    B = read_In_State('State_B.pts');
    targets(:,2) = A(:,1) + (t/t1)*( B(:,1) - A(:,1) );
    targets(:,3) = A(:,2) + (t/t1)*( B(:,2) - A(:,2) );
    
elseif t <= (t1+t2)
    B = read_In_State('State_B.pts');
    C = read_In_State('State_C.pts');
    targets(:,2) = B(:,1) + (t-t1)/(t2) * ( C(:,1) - B(:,1) );
    targets(:,3) = B(:,2) + (t-t1)/(t2) * ( C(:,2) - B(:,2) );
    
elseif t <= (t1+2*t2)
    B = read_In_State('State_B.pts');
    C = read_In_State('State_C.pts');
    targets(:,2) = C(:,1) + (t-t1-t2)/(t2) * ( B(:,1) - C(:,1) );
    targets(:,3) = C(:,2) + (t-t1-t2)/(t2) * ( B(:,2) - C(:,2) );
    
else
    A = read_In_State('State_A.pts');
    B = read_In_State('State_B.pts');
    targets(:,2) = B(:,1) + (t-t1-2*t2)/(t1) * ( A(:,1) - B(:,1) );
    targets(:,3) = B(:,2) + (t-t1-2*t2)/(t1) * ( A(:,2) - B(:,2) );
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
