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
% 	2. Beams (*torsional springs, non-invariant beams*)
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
% FUNCTION: updates the spring attributes!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function springs_info = update_Springs(dt,current_time,xLag,yLag,springs_info)

%springs_info: col 1: starting spring pt (by lag. discretization)
%              col 2: ending spring pt. (by lag. discretization)
%              col 3: spring stiffness
%              col 4: spring resting lengths

%RL = springs_info(:,4); % resting-length vector

freq = 0.875;
%theta=freq*current_time*2*pi;
%phase = sin(theta);

%size(springs_info)

%springs_info(1:3,:)
%pause();

% CHANGE RESTING LENGTH BTWN SIDES OF JELLYFISH BELL
%springs_info(319:end,4) = phase;

springs_info(319:end,4) = abs( sin(0.875*current_time*2*pi) );


% if ( mod(current_time,0.875) < 0.875/2 ) 
%     springs_info(319:end,4) = abs( sin(0.875*current_time*2*pi) );
% else
%     springs_info(319:end,4) = 0;%abs( sin(0.875*current_time*2*pi) );
% end
