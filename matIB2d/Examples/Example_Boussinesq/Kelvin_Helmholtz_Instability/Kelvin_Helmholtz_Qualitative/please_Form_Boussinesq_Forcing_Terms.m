%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  battistn@tcnj.edu
% Date Created: September 9th, 2016
% Institution Created: UNC-CH
% Date Modified: June 26, 2021
% Institution: TCNJ
%
% This code is capable of creating Lagrangian Structures using:
% 	1. Springs
% 	2. Beams (*torsional springs*)
% 	3. Target Points
%	4. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-Tension)")
%
% One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting lengths, etc
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
% 
% If you would like us to add a specific muscle model, please let Nick (battistn@tcnj.edu) know.
%
% If you use this code for the purposes of teaching, research, or recreation please let Nick know as well :)
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Form the Boussinesq forcing matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fBouss_X,fBouss_Y] = please_Form_Boussinesq_Forcing_Terms(exp_Coeff,Nx,Ny,gravity_Info)

    % INPUTS:
    % exp_Coeff:    coefficient of expansion, e.g., thermal expansion
    % Nx,Ny:        grid resolution in x,y directions (NOTE: code not ready for Nx != Ny as of 10/27/2016
    % gravity_Info: 
    %               col 1: flag if considering gravity
    %               col 2: x-component of gravity vector (normalized)
    %               col 3: y-component of gravity vector (normalized)

    g = 9.81;       % gravitational constant
    mat = zeros(Ny,Nx);
        
    %inds = get_HardCode_Inds_Please();
    inds = get_Inds_Please();
    for i=1:length( inds(:,1) )
        xInd = inds(i,1);
        yInd = inds(i,2);
        mat(yInd,xInd) = 1;
    end
    
    fBouss_X = exp_Coeff*g*gravity_Info(2)*mat;
    fBouss_Y = exp_Coeff*g*gravity_Info(3)*mat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in indices for where Boussinesq approximation will take
%           effect in the domain.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function inds = get_Inds_Please()

filename = 'boussinesq.mesh';  %Name of file to read in
fileID = fopen(filename);

% Read in the file, use 'CollectOutput' to gather all similar data together
% and 'CommentStyle' to to end and be able to skip lines in file.
C = textscan(fileID,'%f %f','CollectOutput',1);

fclose(fileID);     %Close the data file.

indices = C{1};    %Stores all read in data in vertices (N+1,2) array

N = length(indices(:,1));  % Total # of indices to affect Boussinesq
inds = zeros(N,2);         % Initialize storage for Lagrangian Pts.

for i=1:N
   inds(i,1) = indices(i,1); %Stores x-values of Lagrangian Mesh
   inds(i,2) = indices(i,2); %Stores y-values of Lagrangian Mesh
  
end
        
