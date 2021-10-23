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
% If you would like us %to add a specific muscle model, please let Nick (nickabattista@gmail.com) know.
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Models the force interaction that is user-defined.
%
%           Note:
%                   1. User has complete control for artibtrary input file
%                      format and number of parameters
%                   2. Those parameters get passed to the matrix called
%                      "general_force"
%                   3. User gets options of including other things into
%                      force function outside of inputted parameters, e.g.,
%                      current_time, dt, previous location of Lag. Pts, etc.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fx_genForce, fy_genForce] = give_Me_General_User_Defined_Force_Densities(ds,Nb,xLag,yLag,xLag_P,yLag_P,dt,current_time,general_force)

    %
    % INPUTS:
    %         ds: Lagragian spacing (defined by ds = 0.5*dx)
    %         Nb: # of Lagrangian Pts.
    %         xLag: current x-Lagrangian coordinate positions
    %         yLag: current y-Lagrangian coordinate positions
    %         xLag_P: previous x-Lagrangian coordinate positions
    %         yLag_P: previous y-Lagrangian coordinate positions
    %         dt: time-step value
    %         current_time: current time in simulation
    %         general_force: matrix containing all data from input file
    %
    
    %
    % NOTE: THIS EXAMPLE CREATES A USER-DEFINED FORCE THAT IS JUST AN
    % ORDINARY LINEAR SPRING
    %
    
    Nsprings = length(general_force(:,1));  % # of Springs
    sp_1 = general_force(:,1);              % Initialize storage for MASTER NODE Spring Connection
    sp_2 = general_force(:,2);              % Initialize storage for SLAVE NODE Spring Connection
    K_Vec = general_force(:,3);             % Stores spring stiffness associated with each spring
    RL_Vec = general_force(:,4);            % Stores spring resting length associated with each spring

    fx = zeros(Nb,1);                 % Initialize storage for x-forces
    fy = fx;                          % Initialize storage for y-forces

    %
    % Loops over all the master-nodes of the springs to compute forces
    %
    for i=1:Nsprings

        id_Master = sp_1(i);          % Master Node index
        id_Slave = sp_2(i);           % Slave Node index
        k_Spring = K_Vec(i);          % Spring stiffness of i-th spring
        L_r = RL_Vec(i);              % Resting length of i-th spring

        dx = xLag(id_Slave) - xLag(id_Master); % x-Distance btwn slave and master node
        dy = yLag(id_Slave) - yLag(id_Master); % y-Distance btwn slave and master node

        sF_x =  k_Spring * ( sqrt( dx^2 + dy^2 ) - L_r ) * ( dx / sqrt(dx^2+dy^2) ); % Compute x-Force
        sF_y =  k_Spring * ( sqrt( dx^2 + dy^2 ) - L_r ) * ( dy / sqrt(dx^2+dy^2) ); % Compute y-Force

        fx(id_Master,1) = fx(id_Master,1) + sF_x;  % Sum total forces for node, i in x-direction (this is MASTER node for this spring)
        fy(id_Master,1) = fy(id_Master,1) + sF_y;  % Sum total forces for node, i in y-direction (this is MASTER node for this spring)

        fx(id_Slave,1) = fx(id_Slave,1) - sF_x;    % Sum total forces for node, i in x-direction (this is SLAVE node for this spring)
        fy(id_Slave,1) = fy(id_Slave,1) - sF_y;    % Sum total forces for node, i in y-direction (this is SLAVE node for this spring)


    end
    
    % Store forces from artibrary force function onto desired Lagrangian Pts.
    fx_genForce = fx; 
    fy_genForce = fy;
