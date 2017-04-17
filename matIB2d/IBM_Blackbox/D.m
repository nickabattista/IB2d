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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Finds CENTERED finite difference approximation to 1ST
% Derivative in specified direction by input, dz, and 'string'. 
% Note: It automatically accounts for periodicity of the domain.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u_z = D(u,dz,string)

% u:      velocity 
% dz:     spatial step in "z"-direction
% string: specifies which 1ST derivative to take (to enforce periodicity)


if strcmp(string,'x')

    len = length(u(1,:));   % LENGTH IN X-DIRECTION
    
    %For periodicity on ends
    u_z(:,1) = ( u(:,2) - u(:,len) ) / (2*dz);
    u_z(:,len)= ( u(:,1) - u(:,len-1) ) / (2*dz);

    %Standard Centered Difference 
    for j=2:len-1
        u_z(:,j) = ( u(:,j+1) - u(:,j-1) ) / (2*dz);
    end

elseif strcmp(string,'y')
    
    len = length(u(:,1));   % LENGTH IN Y-DIRECTION
    
    %For periodicity on ends
    u_z(1,:) = ( u(2,:) - u(len,:) ) / (2*dz);
    u_z(len,:)= ( u(1,:) - u(len-1,:) ) / (2*dz);

    %Standard Centered Difference 
    for j=2:len-1
        u_z(j,:) = ( u(j+1,:) - u(j-1,:) ) / (2*dz);
    end
    
else
    
    fprintf('\n\n\n ERROR IN FUNCTION D FOR COMPUTING 1ST DERIVATIVE\n');
    fprintf('Need to specify which desired derivative, x or y.\n\n\n'); 
       
end


