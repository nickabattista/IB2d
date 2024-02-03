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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Moves Lagrangian Point Positions by doing the integral,
%
%           " xLag_Next = xLag_Prev + dt* int( u(x,t) delta( x - xLag_n ) dX ) "
%
%      NOTE: (i) Lots of old implementation included (but commented out) for teaching purposes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xL_Next, yL_Next] = please_Move_Lagrangian_Point_Positions(mu, u, v, xL_P, yL_P, xL_H, yL_H, x, y, dt, grid_Info,porous_Yes,poroelastic_Yes,poroelastic_info,F_Poro)


% Grid Info
Nx =   grid_Info(1);
Ny =   grid_Info(2);
Lx =   grid_Info(3);
Ly =   grid_Info(4);
dx =   grid_Info(5);
dy =   grid_Info(6);
supp = grid_Info(7);
%Nb =   grid_Info(8);
%ds =   grid_Info(9);


%--------------------------------------------------------------------------------
% Find indices where the delta-function kernels are non-zero for both x and y
%--------------------------------------------------------------------------------
[xInds,yInds] = give_NonZero_Delta_Indices_XY(xL_H, yL_H, Nx, Ny, dx, dy, supp);


%--------------------------------------------------------------------------------
% ReSize the xL_H and yL_H matrices for use in the Dirac-delta function
%        values to find distances between corresponding Eulerian data and them
%--------------------------------------------------------------------------------
xLH_aux = mod(xL_H,Lx); 
yLH_aux = mod(yL_H,Ly); 


%----------------------------------------
%        ORIGINAL IMPLEMENTATION
%  --> Slightly slower
%----------------------------------------
% xL_H_ReSize = []; 
% yL_H_ReSize = [];
% for i=1:supp^2
%    xL_H_ReSize = [xL_H_ReSize xLH_aux];
%    yL_H_ReSize = [yL_H_ReSize yLH_aux];
% end

%------------------------------------------
%   USE REP MAT TO MAKE COPIES QUICKLY
%   --> Slightly faster
%------------------------------------------
xL_H_ReSize = repmat(xLH_aux,1,supp^2);
yL_H_ReSize = repmat(yLH_aux,1,supp^2);


%----------------------------------------------------------------------------
% Checks if only one Lagrangian Point (ensures dimensions line up)
%----------------------------------------------------------------------------
if ( length(xL_P) == 1 )
    % Finds distance between specified Eulerian data and nearby Lagrangian data
    try
        distX = give_Eulerian_Lagrangian_Distance(x(xInds), xL_H_ReSize, Lx);
        distY = give_Eulerian_Lagrangian_Distance(y(yInds)', yL_H_ReSize, Ly);
    catch
        fprintf('\n\n\n - ERROR - \n');
        fprintf('\n\n - ERROR ERROR - \n');
        fprintf('\n\n - ERROR ERROR ERROR - \n');
        fprintf('\n\n - ERROR ERROR ERROR ERROR - \n\n\n');
        error('BLOW UP! (*forces TOO large*) -> try decreasing the time-step or decreasing material property stiffnesses');
    end
else
    
    % Finds distance between specified Eulerian data and nearby Lagrangian data
    try 
        distX = give_Eulerian_Lagrangian_Distance(x(xInds), xL_H_ReSize, Lx);
        distY = give_Eulerian_Lagrangian_Distance(y(yInds), yL_H_ReSize, Ly);
        %size(distX)
    catch
        fprintf('\n\n\n - ERROR - \n');
        fprintf('\n\n - ERROR ERROR - \n');
        fprintf('\n\n - ERROR ERROR ERROR - \n');
        fprintf('\n\n - ERROR ERROR ERROR ERROR - \n\n\n');
        error('BLOW UP! (*forces TOO large*) -> try decreasing the time-step or decreasing material property stiffnesses');
    end
end

%------------------------------------------
% Obtain the Dirac-delta function values
%------------------------------------------
delta_X = give_Delta_Kernel( distX, dx);
delta_Y = give_Delta_Kernel( distY, dy);

%------------------------------------------
% Perform Integral
%------------------------------------------
[move_X, move_Y] = give_Me_Perturbed_Distance(u,v,dx,dy,delta_X,delta_Y,xInds,yInds);

%------------------------------------------
% Update the Lagrangian Point Position
%------------------------------------------
xL_Next = xL_P + (dt) * move_X;
yL_Next = yL_P + (dt) * move_Y;

%----------------------------------------------------------------
% Update the Lagrangian Point Positions with poroelasticity
%----------------------------------------------------------------
if poroelastic_Yes
    %
    % poroelastic_info(:,1): index of poroelastic point
    % poroelastic_info(:,2): Brinkman constant
    %
    xL_Next(poroelastic_info(:,1)) = xL_Next(poroelastic_info(:,1)) + (  1./(mu*poroelastic_info(:,2)) .* F_Poro(:,1) ) * dt;
    yL_Next(poroelastic_info(:,1)) = yL_Next(poroelastic_info(:,1)) + (  1./(mu*poroelastic_info(:,2)) .* F_Poro(:,2) ) * dt;

end

%----------------------------------------------------------------
% TESTING FOR LAG. PTS. MOVING THRU BOUNDARIES
%----------------------------------------------------------------
%xL_NextB = xL_Next;
% if find(xL_NextB>Lx)
%     mat(:,1)=xL_NextB;
%     mat(:,2)=mod(xL_Next,Lx);
%     mat
%     current_time
%     error('problem with lag pts moving through boundary');
% end

%----------------------------------------------------------------
% Shift so that all values are in [0,Lx) or [0,Ly)
%----------------------------------------------------------------
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
        
        % Compute integrand 'stencil' of velocity * delta for each Lagrangian Pt!
        mat_X(i,j) = u(yID,xID)*delta_X(i,j)*delta_Y(i,j);
        mat_Y(i,j) = v(yID,xID)*delta_X(i,j)*delta_Y(i,j);

    end
end

%--------------------------------------------------------------------
% Approximate Integral of Velocity x Delta for each Lagrangian Pt!
%--------------------------------------------------------------------
move_X = sum( mat_X , 2) * (dx*dy);
move_Y = sum( mat_Y , 2) * (dx*dy);
