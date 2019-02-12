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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes a discrete approx. to a 1D Dirac-delta function over a
% specified matrix, x, and spatial step-size, dx. It will have support in
% [x-2dx, x+2dx]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function delta = give_Delta_Kernel(x,dx)

% x:  Values in which the delta function will be evaulated
% dx: Spatial step-size of grid

% Computes Dirac-delta Approximation.
RMAT = abs(x)/dx;

%Initialize delta
delta = RMAT;

%Loops over to find delta approximation
[row,col] = size(x);
for i=1:row
    for j=1:col
        
        r = RMAT(i,j);
        
        % OLD DELTA FUNCTION APPROXIMATION
        %if r<=2
        %    delta(i,j) = 0.25*(1+cos(pi*r/2));
        %else
        %    delta(i,j) = 0;
        %end
        
        % PESKIN 4-PT DISCRETE DELTA FUNCTION
        if r<1
            delta(i,j) = ( (3 - 2*r + sqrt(1 + 4*r - 4*r.*r) ) / (8*dx) );
        elseif ( (r<2) && (r>=1) )
            delta(i,j) = ( (5 - 2*r - sqrt(-7 + 12*r - 4*r.*r) ) / (8*dx) );
        else
            delta(i,j) = 0;
        end
        
    end
end


