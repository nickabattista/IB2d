%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: finds the indices on the Eulerian grid where the 1D Dirac-delta
% kernel is possibly non-zero in BOTH (x,y) directions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%Repeat x-Indices for Non-Zero y-Indices!
xInds = [];
for i=1:supp
   xInds = [xInds xIndsAux]; %Sets up x-INDEX matrix bc we consider BOTH dimensions
end


%Give y-dimension Non-Zero Delta Indices
yIndsAux = give_1D_NonZero_Delta_Indices(yLag, Ny, dy, supp);

%Repeat y-Indices for Non-Zero x-Indices!
yInds = [];
for i=1:supp
    for j=1:supp
        yInds = [yInds yIndsAux(:,i)]; %Sets up y-INDEX matrix bc we consider BOTH dimensions
    end
end





