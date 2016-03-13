%*****************************************************************************%
%***********************************% IB2d %**********************************%
%*****************************************************************************%

% IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
%  	fluid-structure interaction models. This version of the code is based 
% 	off of Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  nick.battista@unc.edu
% Date Created: May 27th, 2015
% Institution: University of North Carolina at Chapel Hill
% Website: http://battista.web.unc.edu
% GitHub: http://www.github.com/nickabattista
% 
% This code is capable of creating Lagrangian Structures using:
%  	1. Springs
%  	2. Beams (*torsional springs)
%  	3. Target Points
% 	4. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-
%                         Tension)")
%         5. Mass Points
% 
% One is able to update those Lagrangian Structure Parameters, e.g., spring constants, resting lengths, etc
%  
% There are a number of built in Examples, mostly used for teaching purposes. 
% 
% If you would like us to add a specific muscle model, please contact Nick (nick.battista@unc.edu) 
% 
% If you use this code for the purposes of teaching, research, or recreation please let Nick know as well :)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: update 'massive' immersed boundary velocity
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mVelocity_h = please_Update_Massive_Boundary_Velocity(dt_step,mass_info,mVelocity,F_Mass_Bnd,gravity_Info)

% dt_step: desired time-step for this position
% mass_info:   col 1: lag index for mass pt
%              col 2: massive x-Lag Value
%              col 3: massive y-Lag Value
%              col 4: 'mass-spring' stiffness parameter
%              col 5: MASS parameter value
% mVelocity    col 1: x-directed Lagrangian velocity
%              col 2: y-directed Lagrangian velocity
% F_Mass_Bnd   col 1: x-directed Lagrangian force
%              col 2: y-directed Lagrangian force
% gravity_Info col 1: flag if considering gravity (0 = NO, 1 = YES)
%              col 2: x-component of gravity vector (NORMALIZED PREVIOUSLY)
%              col 3: y-component of gravity vector (NORMALIZED PREVIOUSLY)

ids = mass_info(:,1);

if gravity_Info(1) == 1
     
    g = 9.80665; %m/s^2
    
    % update x-Velocity
    mVelocity_h(:,1) = mVelocity(:,1) - dt_step * ( F_Mass_Bnd(ids,1)./mass_info(:,5) - g*gravity_Info(2) );

    % update y-Velocity
    mVelocity_h(:,2) = mVelocity(:,2) - dt_step * ( F_Mass_Bnd(ids,2)./mass_info(:,5) - g*gravity_Info(3) );
    
else
    
    % update x-Velocity
    mVelocity_h(:,1) = mVelocity(:,1) - dt_step*F_Mass_Bnd(ids,1)./mass_info(:,5);

    % update y-Velocity
    mVelocity_h(:,2) = mVelocity(:,2) - dt_step*F_Mass_Bnd(ids,2)./mass_info(:,5);
    
end



