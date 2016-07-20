%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled  
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
% If you would like us to add a specific muscle model, please let Nick (nick.battista@unc.edu) know.
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: imports all Eulerian Data at a single step
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,y,Omega,P,uMag,uX,uY,U,V] = import_Eulerian_Data(path,numSim)

% read in Vorticity %
strChoice = 'Omega'; first = 1;
[Omega,x,y] = read_Eulerian_Data_From_vtk(path,numSim,strChoice,first);

% read in Momentum %
strChoice = 'P'; first = 0;
[P,~,~] = read_Eulerian_Data_From_vtk(path,numSim,strChoice,first);

% read in Velocity Magnitude %
strChoice = 'uMag'; first = 0;
[uMag,~,~] = read_Eulerian_Data_From_vtk(path,numSim,strChoice,first);

% read in x-directed Velocity Magnitude %
strChoice = 'uX'; first = 0;
[uX,~,~] = read_Eulerian_Data_From_vtk(path,numSim,strChoice,first);

% read in y-directed Velocity Magnitude %
strChoice = 'uY'; first = 0;
[uY,~,~] = read_Eulerian_Data_From_vtk(path,numSim,strChoice,first);

% read in Velocity Field %
[U,V] = read_Eulerian_Velocity__Field_vtk(path,numSim);

clear strChoice first;

