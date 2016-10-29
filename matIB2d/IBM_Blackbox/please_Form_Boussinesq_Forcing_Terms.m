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
% FUNCTION: update 'massive' immersed boundary position
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
    mat = ones(Ny,Nx);
    fBouss_X = exp_Coeff*g*gravity_Info(2)*mat;
    fBouss_Y = exp_Coeff*g*gravity_Info(3)*mat;
    