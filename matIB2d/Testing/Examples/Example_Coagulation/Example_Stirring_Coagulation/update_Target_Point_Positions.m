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

NP_before = 0;  % # of channel pts before piston targets

%IDs = targets(NP_before+1:end,1);   % Stores Lag-Pt IDs in col vector
xPts= targets(NP_before+1:end,2);    % Original x-Values of x-Target Pts.
%yPts= targets(:,3);                 % Original y-Values of y-Target Pts.
%kStiffs = targets(:,4);             % Stores Target Stiffnesses 

%N_target = length(targets(:,1));    % Gives total number of target pts!

period = 0.05;                      % Time it takes to move to the right
rest = 0.01;
L = 0.5;                             % Distance cylinder moves in each direction
s = L/(period);                      % Speed cylinder moves
t = rem(current_time,2*period+rest); % Time
max_Speed = 15;

%
% Moves wing to the right!
%    
if t <= period

    %offset = tanh(200*current_time)*s*dt;
    
    %coeff = -max_Speed / ( period/4 - period^2/2 );
    %offset = coeff*t*(t-period);
    
    offset = max_Speed*sin( pi*t/period )*dt;
    
    xPts = xPts + offset;
    targets(NP_before+1:end,2) = xPts;

    %for i=1:N_target                    % Loops over all target points!    
    %    xPts(IDs(i)) = xPts(IDs(i)) + tanh(50*current_time)*s*dt;
    %    targets(IDs(i),2) = xPts(IDs(i)); 
    %end
    
elseif t > (period+rest)
    
    %offset = tanh(200*current_time)*s*dt;
    
    %coeff = -max_Speed / ( period/4 - period^2/2 );
    %offset = coeff*(t-period)*(t-2*period);
    
    offset = max_Speed*sin( pi*(t-period-rest)/period )*dt;
    
    xPts = xPts - offset;
    targets(NP_before+1:end,2) = xPts;
    
end

