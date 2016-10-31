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
% One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting %%	lengths, etc
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
% 
% If you would like us %to add a specific muscle model, please let Nick (nick.battista@unc.edu) know.
%
%--------------------------------------------------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: main2d is the function that gets called to run the code. It
%           itself, reads in parameters from the input2d file, and passes 
%           them to the IBM_Driver function to run the simulation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function main2d()

%This is the "main" file, which gets called to run the Immersed Boundary
%Simulation. It reads in all the parameters from "input2d", and sends them
%off to the "IBM_Driver" function to actual perform the simulation


% READ-IN INPUT PARAMETERS %
[params,strName] = give_Me_input2d_Parameters();


% SIMULATION TO RUN %
struct_name = char(strName(1));


% FLUID PARAMETER VALUES STORED %
mu = params(1);      % Dynamic Viscosity
rho = params(2);     % Density

% TEMPORAL INFORMATION VALUES STORED %
T_final = params(3); % Final simulation time
dt = params(4);      % Time-step

% GRID INFO STORED %
grid_Info(1) = params(5);           % # of Eulerian Pts in x-Direction
grid_Info(2) = params(6);           % # of Eulerian Pts in y-Direction 
grid_Info(3) = params(7);           % Length of Eulerian domain in x-Direction
grid_Info(4) = params(8);           % Length of Eulerian domain in y-Direction
grid_Info(5) = params(7)/params(5); % Spatial step-size in x
grid_Info(6) = params(8)/params(6); % Spatial step-size in y
grid_Info(7) = params(9);           % # of pts used in delta-function support (supp/2 in each direction)
grid_Info(8) = params(29);          % Print Dump (How often to plot)
grid_Info(9) = params(30);          % Plot in Matlab? (1=YES,0=NO) 
grid_Info(10) = params(31);         % Plot LAGRANGIAN PTs ONLY in Matlab
grid_Info(11) = params(32);         % Plot LAGRANGIAN PTs + VELOCITY FIELD in Matlab
grid_Info(12) = params(33);         % Plot LAGRANGIAN PTs + VORTICITY colormap in Matlab
grid_Info(13) = params(34);         % Plot LAGRANGIAN PTs + MAGNITUDE OF VELOCITY colormap in Matlab
grid_Info(14) = params(35);         % Plot LAGRANGIAN PTs + PRESSURE colormap in Matlab


% MODEL STRUCTURE DATA STORED %
model_Info(1) = params(10);         % Springs: 0 (for no) or 1 (for yes) 
model_Info(2) = params(11);         % Update_Springs: 0 (for no) or 1 (for yes)
model_Info(3) = params(12);         % Target_Pts: 0 (for no) or 1 (for yes)
model_Info(4) = params(13);         % Update_Target_Pts: 0 (for no) or 1 (for yes)
model_Info(5) = params(14);         % Beams: 0 (for no) or 1 (for yes)
model_Info(6) = params(15);         % Update_Beams: 0 (for no) or 1 (for yes)
model_Info(7) = params(16);         % Muscle Activation (Length/Tension-Hill Model): 0 (for no) or 1 (for yes)
model_Info(8) = params(17);         % Muscle Activation 3-ELEMENT HILL MODEL w/ Length-Tension/Force-Velocity: 0 (for no) or 1 (for yes)
model_Info(9) = params(18);         % Arbirtary External Force Onto Fluid Grid: 0 (for no) or 1 (for yes)
model_Info(10) = params(19);        % Tracer Particles: 0 (for no) or 1 (for yes)
model_Info(11)= params(20);         % Mass Points: 0 (for no) or 1 (for yes)
model_Info(12)= params(21);         % Gravity: 0 (for no) or 1 (for yes)
model_Info(13)= params(22);         % x-Component of Gravity vector
model_Info(14)= params(23);         % y-Component of Gravity Vector
model_Info(15)= params(24);         % Porous Media: 0 (for no) or 1 (for yes)
model_Info(16)= params(25);         % Background Concentration Gradient: 0 (for no) or 1 (for yes)
model_Info(17)= params(26);         % Electrophysiology Model (FitzHugh-Nagumo)
model_Info(18)= params(27);         % Damped Springs: 0 (for no) or 1 (for yes)
model_Info(19)= params(28);         % Update_Damped_Springs: 0 (for no) or 1 (for yes)

% Path Reference to where Driving code is found %
addpath('../../../IBM_Blackbox/');

%-%-%-% DO THE IMMERSED BOUNDARY SOLVE!!!!!!!! %-%-%-%
[X, Y, U, V, xLags, yLags] = IBM_Driver(struct_name, mu, rho, grid_Info, dt, T_final, model_Info);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to read in input files from "input2d" -> renames all
% quantities appropriately, just as they are in input file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [params,struct_name] = give_Me_input2d_Parameters()

filename= 'input2d';  %Name of file to read in

fileID = fopen(filename); %Opens file for 'textscan' function
    % Read in the file, use 'CollectOutput' to gather all similar data together
    % and 'CommentStyle' to to end and be able to skip lines in file.
    C1 = textscan(fileID,'%s %s %f','CollectOutput',1,'CommentStyle','%');
    C2 = textscan(fileID,'%s %s %s','CollectOutput',1,'CommentStyle','%');
fclose(fileID);

params = C1{2};

struct_name = C2{1,1};

