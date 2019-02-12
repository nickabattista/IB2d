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
% FUNCTION finds the indices on the Eulerian grid where the 1D Dirac-delta
% kernel is possibly non-zero is x-dimension.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function indices = give_1D_NonZero_Delta_Indices(lagPts_j, N, dx, supp)

%lagPts_j: matrix of lagrangian pts for specific coordinate, j= x or y.
%N:        # spatial resolution of Eulerian grid in each dimension
%dx:       Spatial step-size on Eulerian (fluid) grid
%supp:     Size of support of the Dirac-delta kernel (should be even)


% Finds the index of the lower left Eulerian pt. to Lagrangian pt..
ind_Aux = floor(lagPts_j/dx + 1);

% Get all the different x indices that must be considered.
indices = [];
for i=1:supp
    indices = [indices ind_Aux]; 
end
%
for i=1:supp
    indices(:,i) = indices(:,i) + -supp/2+1+(i-1); 
end

% Translate indices between {1,2,..,N}
indices = mod(indices-1,N) + 1;
