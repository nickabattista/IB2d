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

function [fX_Lag,fY_Lag,fLagMag,fLagNorm,fLagTan] = import_Lagrangian_Force_Data(path,numSim)

% read in Mag. of Force %
strChoice = 'fX_Lag'; 
fX_Lag =  read_Force_Scalar_Data_From_vtk(path,numSim,strChoice);

% read in Mag. of Tangential Force %
strChoice = 'fY_Lag'; 
fY_Lag = read_Force_Scalar_Data_From_vtk(path,numSim,strChoice);

% read in Mag. of Force %
strChoice = 'fMag'; 
fLagMag =  read_Force_Scalar_Data_From_vtk(path,numSim,strChoice);

% read in Mag. of Tangential Force %
strChoice = 'fTan'; 
fLagNorm = read_Force_Scalar_Data_From_vtk(path,numSim,strChoice);

% read in Mag. of Normal Force %
strChoice = 'fNorm';
fLagTan =  read_Force_Scalar_Data_From_vtk(path,numSim,strChoice);
 

clear strChoice first;

