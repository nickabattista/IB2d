%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  nickabattista@gmail.com
% Date Created: May 27th, 2015
% Institution: UNC-CH
%
% This code is capable of creating Lagrangian Structures using:
% 	1. Springs
% 	2. Beams (*torsional springs)
% 	3. Target Points
%	4. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-Tension)")
%
% One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting lengths, etc
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
% 
% If you would like us %to add a specific muscle model, please let Nick (nickabattista@gmail.com) know.
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Computes the source or sink term for concentration at the moving boundary
%          
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Fs] = please_Find_Source_For_Concentration(dt, current_time, xLag, yLag, x, y, grid_Info, model_Info, k, C, ci, flag_Geo_Connect, geo_Connect_MAT)

%
% dt:             time step
% current_time:   current time in simulation
% xLag:           x positions of Lagrangian structure
% yLag:           y positions of Lagrangian structure  
% x:              x positions on Eulerian grid
% y:              y positions on Eulerian grid
% grid_Info:      holds lots of geometric pieces about grid / simulations
% model_Info:     1 if constant model (int k delta(x-xLag)delta*(y-yLag)ds) or 2 if nonconstant
%                 model (int k(ci-C) delta(x-xLag)delta*(y-yLag)ds)  
% k:              Rate of concentration creation
% C:              Concentration 
% ci:             Concentration saturation limit
% c_Per:          Concentration periodic boundary 
       
% Grid Info %
Nx =    grid_Info(1); % # of Eulerian pts. in x-direction
Ny =    grid_Info(2); % # of Eulerian pts. in y-direction
Lx =    grid_Info(3); % Length of Eulerian grid in x-coordinate
Ly =    grid_Info(4); % Length of Eulerian grid in y-coordinate
dx =    grid_Info(5); % Spatial-size in x
dy =    grid_Info(6); % Spatial-size in y
supp =  grid_Info(7); % Delta-function support
Nb =    grid_Info(8); % # of Lagrangian pts. 
ds =    grid_Info(9); % Lagrangian spacing


% ReSize the xL_H and yL_H matrices for use in the Dirac-delta function
%        values to find distances between corresponding Eulerian data and them
xLH_aux = mod(xLag,Lx); xL_H_ReSize = [];
yLH_aux = mod(yLag,Ly); yL_H_ReSize = [];
for i=1:supp^2
   xL_H_ReSize = [xL_H_ReSize xLH_aux];
   yL_H_ReSize = [yL_H_ReSize yLH_aux];
end


% Find indices where the delta-function kernels are non-zero for both x and y.
[xInds,yInds] = give_NonZero_Delta_Indices_XY(xLag, yLag, Nx, Ny, dx, dy, supp);



if ( length(xLag) == 1 )
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

    catch
        fprintf('\n\n\n - ERROR - \n');
        fprintf('\n\n - ERROR ERROR - \n');
        fprintf('\n\n - ERROR ERROR ERROR - \n');
        fprintf('\n\n - ERROR ERROR ERROR ERROR - \n\n\n');
        error('BLOW UP! (*forces TOO large*) -> try decreasing the time-step or decreasing material property stiffnesses');
    end
end



% Obtain the Dirac-delta function values.
delta_X = give_Delta_Kernel( distX, dx);
delta_Y = give_Delta_Kernel( distY, dy);

% Perform Integral finding concentration on the Lagrangian boundary
CL = give_Lagrangian_Concentration(C,dx,dy,delta_X,delta_Y,xInds,yInds);

    if model_Info==1 
    fs = zeros(Nb,1)+k;        %constant
    elseif model_Info==2

% Perform Integral
%CL = give_Lagrangian_Concentration(C,dx,dy,delta_X,delta_Y,xInds,yInds);

    fs=k*(ci-CL);         % limited

    elseif model_Info ==3

    fs=k*CL;
    else
  
    'Invalid source model!' 
    end 
    
% Give me delta-function approximations!
[delta_X, delta_Y] = give_Me_Delta_Function_Approximations_For_Force_Calc(x,y,grid_Info,xLag,yLag);

%------------------------------------------------------------
% Transform the force density vectors into diagonal matrices
%------------------------------------------------------------
fsds = zeros(Nb,Nb);

if flag_Geo_Connect % Compute -actual- distances between 
                    % neighboring Lagrangian pts in geometry
                    
    for i=1:Nb
        
        % compute real distances btwn LAG_i and Attached Points
        %                           and sum distances together
        
        % find all Lag Pts that are Lag_i's geometric neighbors
        indsVec = find(geo_Connect_MAT(:,1)==i);
    
        % loop to (1) distances between Lag_i and each neighbor
        %         (2) sum all such distances together
        ds_sum = 0;
        for j=1:length(indsVec)
            id1 = geo_Connect_MAT( indsVec(j),1);
            id2 = geo_Connect_MAT( indsVec(j),2);
            ds_sum = ds_sum + sqrt( ( xLag(id1,1)-xLag(id2,1) )^2 + ( yLag(id1,1)-yLag(id2,1) )^2 );     
        end
        
        % Note: (1) 0.5 coming from Trapezoid Rule
        %       (2) ds_sum is sum of distances to neighboring pts
        fsds(i,i) = 0.5*fs(i,1)*ds_sum; 
    end
        
else % Use Peskin's Constant 'ds' assumption
    for i=1:Nb
       fsds(i,i) = fs(i,1)*ds; 
    end
end

% Find source term on grids by approximating the line integral, 
%       F(x,y) = int{ CL(s) delta(x - xLag(s)) delta(y - yLag(s)) ds }
Fs = delta_Y * fsds * delta_X;

% debugging
%save(sprintf('fs_%f.mat',current_time),'Fs')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION computes the Delta-Function Approximations 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [delta_X delta_Y] = give_Me_Delta_Function_Approximations_For_Force_Calc(x,y,grid_Info,xLag,yLag)

% Grid Info
Nx =   grid_Info(1);
Ny =   grid_Info(2);
Lx =   grid_Info(3);
Ly =   grid_Info(4);
dx =   grid_Info(5);
dy =   grid_Info(6);
supp = grid_Info(7);
Nb =   grid_Info(8);

% Find the indices of the points (xi, yj) where the 1D delta functions are non-zero in EULERIAN FRAME
indX = give_1D_NonZero_Delta_Indices(xLag, Nx, dx, supp);
indY = give_1D_NonZero_Delta_Indices(yLag, Ny, dy, supp)';

% Matrix of possible indices, augmented by "supp"-copies to perform subtractions later in LAGRANGIAN FRAME
indLagAux = [1:1:Nb]';
ind_Lag = [];
for i=1:supp
   ind_Lag = [ind_Lag indLagAux]; 
end


% Compute distance between Eulerian Pts. and Lagrangian Pts. by passing correct indices for each
try
    distX = give_Eulerian_Lagrangian_Distance(x(indX),xLag(ind_Lag),Lx);
    distY = give_Eulerian_Lagrangian_Distance(y(indY),yLag(ind_Lag'),Ly);
catch
    fprintf('\n\n\n - ERROR - \n');
    fprintf('\n\n - ERROR ERROR - \n');
    fprintf('\n\n - ERROR ERROR ERROR - \n');
    fprintf('\n\n - ERROR ERROR ERROR ERROR - \n\n\n');
    error('BLOW UP! (*SOURCE TOO large*) -> try decreasing the time-step or diffusion coefficient');
end

% Initialize delta_X and delta_Y matrices for storing delta-function info for each Lag. Pt.
delta_X = zeros(Nb, Nx);
delta_Y = zeros(Ny, Nb);

delta_1D_x = give_Delta_Kernel(distX, dx);
delta_1D_y = give_Delta_Kernel(distY, dy);


[row,col] = size(ind_Lag);
for i=1:row
    for j=1:col
        
        % Get Eulerian/Lagrangian indices to use for saving non-zero delta-function values
        xID = indX(i,j);
        indy= ind_Lag(i,j);
        yID = indY(j,i);
        
        % Store non-zero delta-function values into delta_X / delta_Y matrices at correct indices!
        delta_X(indy,xID) = delta_1D_x(i,j);
        delta_Y(yID,indy) = delta_1D_y(j,i);
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Computes the integral to find concentration on the Lagrangian points
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C_Lag] = give_Lagrangian_Concentration(c,dx,dy,delta_X,delta_Y,xInds,yInds)

% c:        x-component of velocity
% delta_X:  values of Dirac-delta function in x-direction
% delta_Y:  values of Dirac-delta function in y-direction
% xInds:    x-Indices on fluid grid
% yInds:    y-Indices on fluid grid


[row,col] = size(xInds);
mat = zeros(row,col);  % Initialize matrix for storage
for i=1:row
    for j=1:col
        
        % Get Eulerian indices to use for concentration grid 
        xID = xInds(i,j);
        yID = yInds(i,j);
        
        % Compute integrand 'stencil' of concentration delta for each Lagrangian Pt!
        mat(i,j) = c(yID,xID)*delta_X(i,j)*delta_Y(i,j);

    end
end

% Approximate Integral of concentration x Delta for each Lagrangian Pt!
C_Lag = sum( mat, 2) * (dx*dy);

