%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  nickabattista@gmail.com
% Date Created: October 8th, 2018
% Institution: UNC-CH
%
% This code is capable of creating Lagrangian Structures using:
% 	1. Springs
% 	2. Beams (*torsional springs / non-invariant beams *)
% 	3. Target Points
%	4. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-Tension)")
%
% One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting lengths, etc
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
% 
% If you would like us %to add a specific muscle model, please let Nick (nickabattista@gmail.com) know.
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: helps input previous .vtk information for restarting a
%           simulation that has ended because of power failure, etc.
%
%      NOTE: for restart protocol, need to have .vtk data for:
%                       1. lagPts (Lagrangian positions)
%                       2. u (velocity field)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [current_time,cter,ctsave,U,V,xLag,yLag,xLag_P,yLag_P,path_to_data] = help_Me_Restart(dt)

% 
% NEEDS TO BE HARDCODED PER SIMULATION BASIS
% 
ctsave = 3;              % Last time-step of data saved (# at end of .vtk file);
print_dump = 20;         % Print_dump interval as given in input2d

% Path to Simulation Data (e.g., including viz_IB2d folder)
path_to_data = '/Users/battistn/Desktop/IB2d/matIB2d/Examples/Example_Test_Restart_Protocol/Poroelastic_Rubberband_Restart/viz_IB2d/';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------- DO NOT CHANGE BELOW --------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[current_time,cter,ctsave,U,V,xLag,yLag,xLag_P,yLag_P] = pass_Back_Data_For_Restart(dt,ctsave,print_dump,path_to_data);


