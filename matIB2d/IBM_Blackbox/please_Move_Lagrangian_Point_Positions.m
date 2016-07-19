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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Moves Lagrangian Point Positions by doing the integral,
%
%           " xLag_Next = xLag_Prev + dt* int( u(x,t) delta( x - xLag_n ) dX ) "
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xL_Next, yL_Next] = please_Move_Lagrangian_Point_Positions(u, v, xL_P, yL_P, xL_H, yL_H, x, y, dt, grid_Info,porous_Yes)


% Grid Info
Nx =   grid_Info(1);
Ny =   grid_Info(2);
Lx =   grid_Info(3);
Ly =   grid_Info(4);
dx =   grid_Info(5);
dy =   grid_Info(6);
supp = grid_Info(7);
Nb =   grid_Info(8);
ds =   grid_Info(9);


% Find indices where the delta-function kernels are non-zero for both x and y.
[xInds,yInds] = give_NonZero_Delta_Indices_XY(xL_H, yL_H, Nx, Ny, dx, dy, supp);

% ReSize the xL_H and yL_H matrices for use in the Dirac-delta function
%        values to find distances between corresponding Eulerian data and them
xLH_aux = mod(xL_H,Lx); xL_H_ReSize = [];
yLH_aux = mod(yL_H,Ly); yL_H_ReSize = [];
for i=1:supp^2
   xL_H_ReSize = [xL_H_ReSize xLH_aux];
   yL_H_ReSize = [yL_H_ReSize yLH_aux];
end

% Finds distance between specified Eulerian data and nearby Lagrangian data
distX = give_Eulerian_Lagrangian_Distance(x(xInds), xL_H_ReSize, Lx);
distY = give_Eulerian_Lagrangian_Distance(y(yInds), yL_H_ReSize, Ly);

% Obtain the Dirac-delta function values.
delta_X = give_Delta_Kernel( distX, dx);
delta_Y = give_Delta_Kernel( distY, dy);

% Perform Integral
[move_X, move_Y] = give_Me_Perturbed_Distance(u,v,dx,dy,delta_X,delta_Y,xInds,yInds);

% Update the Lagrangian Point Position.
xL_Next = xL_P + (dt) * move_X;
yL_Next = yL_P + (dt) * move_Y;


% Shift so that all values are in [0,Lx or Ly).
if porous_Yes == 0
    xL_Next = mod(xL_Next, Lx);
    yL_Next = mod(yL_Next, Ly);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Computes the integral to move each Lagrangian Pt!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [move_X, move_Y] = give_Me_Perturbed_Distance(u,v,dx,dy,delta_X,delta_Y,xInds,yInds)

% u:        x-component of velocity
% v:        y-component of velocity
% delta_X:  values of Dirac-delta function in x-direction
% delta_Y:  values of Dirac-delta function in y-direction
% xInds:    x-Indices on fluid grid
% yInds:    y-Indices on fluid grid


[row,col] = size(xInds);
mat_X = zeros(row,col);  % Initialize matrix for storage
mat_Y = zeros(row,col);  % Initialize matrix for storage
for i=1:row
    for j=1:col
        
        % Get Eulerian indices to use for velocity grids, u and 
        xID = xInds(i,j);
        yID = yInds(i,j);
        
        % Compute integrand 'stencil' of velocity x delta for each Lagrangian Pt!
        mat_X(i,j) = u(yID,xID)*delta_X(i,j)*delta_Y(i,j);
        mat_Y(i,j) = v(yID,xID)*delta_X(i,j)*delta_Y(i,j);

    end
end

% Approximate Integral of Velocity x Delta for each Lagrangian Pt!
move_X = sum( mat_X , 2) * (dx*dy);
move_Y = sum( mat_Y , 2) * (dx*dy);
