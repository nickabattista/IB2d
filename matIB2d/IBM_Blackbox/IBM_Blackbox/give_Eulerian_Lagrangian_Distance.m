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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION distance between Eulerian grid data, x, and Lagrangian grid data, y, 
%          at specifed pts typically and makes sure the distance are [0,L] accordingly.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function distance = give_Eulerian_Lagrangian_Distance(x, y, L)

%x,y: two matrices that you find the distance between (x-typically Eulerian data, y-typically Lagrangian data)
%L:   length of domain, i.e., [0,L]

%[row,col] = size(x);
distance = abs( x - y );

% VECTORIZED function calculations for speedup:
distance = min( distance, L-distance );

% SLOWER (for-loop, non-vectorized computation)
% for i=1:row
%     for j=1:col
%         distance(i,j) = min( distance(i,j), L-distance(i,j) ); %Note: need to make sure that taking correct value
%     end
% end

