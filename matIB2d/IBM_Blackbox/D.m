%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  battistn[@]tcnj[.]edu
% 
% IB2d was Created: May 27th, 2015 at UNC-CH
%
% This code is capable of creating Lagrangian Structures using:
% 	1. Springs
% 	2. Beams (*torsional springs)
% 	3. Target Points
% 	4. Muscle-Model (combined Force-Length-Velocity model, "Hill+(Length-Tension)")
%   .
%   .
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Finds CENTERED finite difference approximation to 1ST
%           Derivative in specified direction by input, dz, and 'string'. 
%
%   Note: (1) It automatically accounts for periodicity of the domain.
%         (2) Old implementation included (but commented out) for teaching purposes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u_z = D(u,dz,string)

% u:      velocity 
% dz:     spatial step in "z"-direction
% string: specifies which 1ST derivative to take (to enforce periodicity)

% Initialize Storage
u_z = zeros(size(u));

if strcmp(string,'x')

    len = length(u(1,:));   % LENGTH IN X-DIRECTION
    
    %--------------------------------------------    
    % ORIGINAL Non-vectorized version (slower)   
    %--------------------------------------------
    %     %For periodicity on ends
    %     u_z(:,1) = ( u(:,2) - u(:,len) ) / (2*dz);
    %     u_z(:,len)= ( u(:,1) - u(:,len-1) ) / (2*dz);
    % 
    %     %Standard Centered Difference 
    %     for j=2:len-1
    %         u_z(:,j) = ( u(:,j+1) - u(:,j-1) ) / (2*dz);
    %     end

    % For periodicity on ends
    u_z(:,1) = ( u(:,2) - u(:,len) ) / (2*dz);
    
    % Central Differences (interior of domain)
    u_z(:,2:len-1) = ( u(:,3:len) - u(:,1:len-2) ) / (2*dz);
    
    % For periodicity on ends
    u_z(:,len)= ( u(:,1) - u(:,len-1) ) / (2*dz);


elseif strcmp(string,'y')
    
    len = length(u(:,1));   % LENGTH IN Y-DIRECTION

    %--------------------------------------------    
    % ORIGINAL Non-vectorized version (slower)   
    %--------------------------------------------    
    %     %For periodicity on ends
    %     u_z(1,:) = ( u(2,:) - u(len,:) ) / (2*dz);
    %     u_z(len,:)= ( u(1,:) - u(len-1,:) ) / (2*dz);
    % 
    %     %Standard Centered Difference 
    %     for j=2:len-1
    %         u_z(j,:) = ( u(j+1,:) - u(j-1,:) ) / (2*dz);
    %     end

    % For periodicity on ends
    u_z(1,:) = ( u(2,:) - u(len,:) ) / (2*dz);
    
    % Central Differences (interior of domain)
    u_z(2:len-1,:) = ( u(3:len,:) - u(1:len-2,:) ) / (2*dz);
    
    % For periodicity on ends
    u_z(len,:)= ( u(1,:) - u(len-1,:) ) / (2*dz);

else
    
    fprintf('\n\n\n ERROR IN FUNCTION D FOR COMPUTING 1ST DERIVATIVE\n');
    fprintf('Need to specify which desired derivative, x or y.\n\n\n'); 
       
end


