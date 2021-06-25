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
% If you would like us %to add a specific muscle model, please let Nick (nick.battista@unc.edu) know.
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: updates the resting lengths and/or spring stiffness
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function springs = update_Springs(dt,current_time,xLag,yLag,springs)

% dt: time-step
% current_time: current time of simulation
% xLag: Vector of all xLagrangian Pts        
% springs: col 1: starting spring pt (by lag. discretization)
%          col 2: ending spring pt. (by lag. discretization)
%          col 3: spring stiffness
%          col 4: spring resting lengths

%s_1 = springs(:,1);    % LEADER Node Indices
%s_2 = springs(:,2);    % FOLLOWER Node Indices
k_Vec = springs(:,3);  % Springs Stiffnesses Vector
%RL_Vec= springs(:,4);  % Resting Lengths Vectors

%N = length(xLag);            % Gives total number of Lagrangian pts!
Ns = length( springs(:,1) ); % Gives total number of springs!
Ly = 1.0;
Ny = 128;
ds = 0.5*Ly/Ny;

distVec = abs( sqrt( ( xLag(2:Ns+1)-xLag(1:Ns) ).^2 + ( yLag(2:Ns+1)-yLag(1:Ns) ).^2 ) );

for i=1:Ns
    if distVec(i) > 5*ds
        k_Vec(i) = 0;
        %fprintf('\nSpring %d Connection Broken at Time = %d\n',i,current_time);
    end
end

springs(:,3) = k_Vec; % Update the springs_info

