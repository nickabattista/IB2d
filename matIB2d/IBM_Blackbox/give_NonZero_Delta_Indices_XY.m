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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: finds the indices on the Eulerian grid where the 1D Dirac-delta
%           kernel is possibly non-zero in BOTH (x,y) directions
%
%      NOTE: (i) Lots of old implementation included (but commented out) for teaching purposes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xInds,yInds] = give_NonZero_Delta_Indices_XY(xLag, yLag, Nx, Ny, dx, dy, supp)

%xLag: gives x-coordinate of Lagrangian position
%yLag: gives y-coordinate of Lagrangian position
%Nx:   # of Eulerian grid pts. in x-dimension
%Ny:   # of Eulerian grid pts. in y-dimension
%dx:   spatial-step along x-dimension of Eulerian grid
%dy:   spatial-step along y-dimension of Eulerian grid
%supp: size of support of the Dirac-delta kernel (should be even)


%Give x-dimension Non-Zero Delta Indices
xIndsAux = give_1D_NonZero_Delta_Indices(xLag, Nx, dx, supp);

%---------------------------------------------
%        ORIGINAL IMPLEMENTATION
% Repeat x-Indices for Non-Zero y-Indices!
%---------------------------------------------
% xInds = [];
% for i=1:supp
%    xInds = [xInds xIndsAux]; %Sets up x-INDEX matrix bc we consider BOTH dimensions
% end

%---------------------------------------------
%      TRY TO SPEED UP IMPLEMENTATION
% Repeat x-Indices for Non-Zero y-Indices!
% (note: not a real substantial CPU savings..)
%---------------------------------------------
xInds = repmat(xIndsAux,1,supp);

%---------------------------------------------
% Give y-dimension Non-Zero Delta Indices
%---------------------------------------------
yIndsAux = give_1D_NonZero_Delta_Indices(yLag, Ny, dy, supp);

%---------------------------------------------
%        ORIGINAL IMPLEMENTATION
% Repeat y-Indices for Non-Zero x-Indices!
%---------------------------------------------
% yInds = [];
% for i=1:supp
%     for j=1:supp
%         yInds = [yInds yIndsAux(:,i)]; %Sets up y-INDEX matrix bc we consider BOTH dimensions
%     end
% end

%---------------------------------------------
%     SLIGHTLY FASTER IMPLEMENTATION
% Repeat y-Indices for Non-Zero y-Indices!
%---------------------------------------------
yInds = zeros(length(yIndsAux(:,1)),supp^2);
for i=1:supp
    id1 = 1+(i-1)*supp;
    id2 = i*supp;
    yInds(:,id1:id2) =  repmat(yIndsAux(:,i),1,supp);
end


