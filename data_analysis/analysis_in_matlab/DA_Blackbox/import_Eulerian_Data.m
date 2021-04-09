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

function [x,y,Omega,P,uMag,uX,uY,U,V,Fx,Fy,C] = import_Eulerian_Data_2(path,numSim,Eulerian_Flags)

    %
    % EULERIAN FLAGS FOR WHAT GETS SPIT OUT %
    % 
    %Eulerian_Flags(1):   OMEGA
    %Eulerian_Flags(2):   PRESSURE
    %Eulerian_Flags(3):   uMAG
    %Eulerian_Flags(4):   uX (mag. x-component of velocity)
    %Eulerian_Flags(5):   uY (mag. x-component of velocity)
    %Eulerian_Flags(6):   uVEC (x,y-components of velocity: U,V)
    %Eulerian_Flags(7):   Fx (x-component of force )
    %Eulerian_Flags(8):   Fy (y-component of force)
    %

    
analysis_path = pwd;
    
% read in Vorticity %
if Eulerian_Flags(1)
    strChoice = 'Omega'; first = 1;
    [Omega,x,y] = read_Eulerian_Data_From_vtk(path,numSim,strChoice,first);
else
    Omega=[];
end


% read in Pressure %
if Eulerian_Flags(2)
    strChoice = 'P'; first = 1;
    [P,x,y] = read_Eulerian_Data_From_vtk(path,numSim,strChoice,first);
else
    P=[];
end


% read in Velocity Magnitude %
if Eulerian_Flags(3)
    strChoice = 'uMag'; first = 1;
    [uMag,x,y] = read_Eulerian_Data_From_vtk(path,numSim,strChoice,first);
else
    uMag=[];
end

% read in x-directed Velocity Magnitude %
if Eulerian_Flags(4)
    strChoice = 'uX'; first = 1;
    [uX,x,y] = read_Eulerian_Data_From_vtk(path,numSim,strChoice,first);
else
    uX=[];
end

% read in y-directed Velocity Magnitude %
if Eulerian_Flags(5)
    strChoice = 'uY'; first = 1;
    [uY,x,y] = read_Eulerian_Data_From_vtk(path,numSim,strChoice,first);
else
    uY=[];
end

% read in x-directed Forces %
if Eulerian_Flags(7)
    strChoice = 'Fx'; first = 1;
    [Fx,x,y] = read_Eulerian_Data_From_vtk(path,numSim,strChoice,first);
else
    Fx=[];
end

% read in y-directed Forces %
if Eulerian_Flags(8)
    strChoice = 'Fy'; first = 1;
    [Fy,x,y] = read_Eulerian_Data_From_vtk(path,numSim,strChoice,first);
else
    Fy=[];
end

% read in Concentration %
if Eulerian_Flags(9)
    strChoice = 'concentration'; first = 1;
    [C,x,y] = read_Eulerian_Data_From_vtk(path,numSim,strChoice,first);
else
    C=[];
end

% read in Velocity Field %
if Eulerian_Flags(6)
    [U,V] = read_Eulerian_Velocity_Field_vtk(path,numSim);
else
    U=[];
    V=[];
end

% Default for x,y values
if max(Eulerian_Flags([1:5,7:9]))==0
    x=[];
    y=[];
end

cd(analysis_path);

clear analysis_path;

clear strChoice first;

