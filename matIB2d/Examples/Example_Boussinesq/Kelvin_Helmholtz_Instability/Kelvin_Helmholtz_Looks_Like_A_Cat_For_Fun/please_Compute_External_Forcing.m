%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  battistn@tcnj.edu
% Date Created: September 9th, 2016
% Institution Created: UNC-CH
% Date Modified: June 26, 2021
% Institution: TCNJ
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
% If you would like us to add a specific muscle model, please let Nick (battistn@tcnj.edu) know.
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Computes the components of the force term in Navier-Stokes from
%           arbitrary external forces, i.e., external force to get desired
%           velocity profile on fluid grid
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Fx, Fy, first, inds] = please_Compute_External_Forcing(dt,current_time,x,y, grid_Info, uX, uY, first, inds)

%
% dt:           time-step 
% current_time: Current time of simulation (in seconds)
% x:            x-Eulerian pts
% y:            y-Eulerian pts
% grid_Info:    holds lots of geometric pieces about grid / simulations
% uX:           x-Velocity on Eulerian Grid
% uY:           y-Velocity on Eulerian Grid


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


% Compute Where You Want to Apply Force
%
xMin1 = 0.0;
xMax1 = Lx;
%
% xMin2 = 0.23;
% xMax2 = 0.27;
% %
% xMin3 = 0.48;
% xMax3 = 0.52;
% %
% xMin4 = 0.73;
% xMax4 = 0.77;
% %
% xMin5 = 0.98;
% xMax5 = 1.0;
%
yMin = 0;
yMax = Ly;

% Stiffness for Arbitrary External Force to Fluid Grid
kStiff = 1e4;

% Width of Channel
w = 0.3;

% Max Velocity Desired
uMax = 1;

if first == 1
    % IF NO SECTIONS
    inds = give_Me_Indices_To_Apply_Force(x,y,xMin1,xMax1,yMin,yMax);
    
    % IF SECTIONS
    %inds1 = give_Me_Indices_To_Apply_Force(x,y,xMin1,xMax1,yMin,yMax);
    %inds2 = give_Me_Indices_To_Apply_Force(x,y,xMin2,xMax2,yMin,yMax);
    %inds3 = give_Me_Indices_To_Apply_Force(x,y,xMin3,xMax3,yMin,yMax);
    %inds4 = give_Me_Indices_To_Apply_Force(x,y,xMin4,xMax4,yMin,yMax);
    %inds5 = give_Me_Indices_To_Apply_Force(x,y,xMin5,xMax5,yMin,yMax);
    %inds = [inds1; inds2; inds3; inds4; inds5];
    first = 0;
end

% Compute External Forces from Desired Target Velocity
[fx, fy] = give_Me_Velocity_Target_External_Force_Density(current_time,dx,dy,x,y,Nx,Ny,Lx,Ly,uX,uY,kStiff,w,uMax,inds);
    
% Compute Total External Forces
Fx = fx;
Fy = fy;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes indices for exerting forces in specified places on fluid grid 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function inds = give_Me_Indices_To_Apply_Force(x,y,xMin,xMax,yMin,yMax)

j=1; noMinYet = 1;
while noMinYet
    
    if ( x(j) >= xMin )
        iX_min = j;
        noMinYet = 0;
    end
    j=j+1;
end

j=length(x); noMaxYet = 1;
while noMaxYet
    
    if ( x(j) <= xMax )
        iX_max = j;
        noMaxYet = 0;
    end
    j=j-1;
end

j=1; noMinYet = 1;
while noMinYet
    
    if ( y(j) >= yMin )
        iY_min = j;
        noMinYet = 0;
    end
    j=j+1;
end

j=length(y); noMaxYet = 1;
while noMaxYet
    
    if ( y(j) <= yMax )
        iY_max = j;
        noMaxYet = 0;
    end
    j=j-1;
end

iX_Vec = iX_min:1:iX_max;
iY_Vec = iY_min:1:iY_max;

n = 1;
for i=1:length(iX_Vec)
    for j=1:length(iY_Vec)
        inds(n,1) = iX_Vec(i);
        inds(n,2) = iY_Vec(j);
        n = n+1; 
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes the External Force Densities! 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fx_exts, fy_exts] = give_Me_Velocity_Target_External_Force_Density(t,dx,dy,x,y,Nx,Ny,Lx,Ly,uX,uY,kStiff,w,Umax,inds)

% t:  current time in simulation
% Nx: # of nodes in x-direction on Eulerian grid
% Ny: # of nodes in y-direction on Eulerian grid
% uX: x-Velocity on Eulerian grid
% uY: y-Velocity on Eulerian grid
% kStiff: stiffness parameter
% inds: indices on the fluid grid for where to apply the arbitrary external force


fx = zeros(Ny,Nx);         % Initialize storage for x-force density from EXTERNAL FORCES
fy = fx;                   % Initialize storage for y-force density from EXTERNAL FORCES

%if t<0.01
    for n=1:length(inds(:,1))
        
        i = inds(n,1);
        j = inds(n,2);

        [uX_Tar] = please_Give_Target_Velocity(t,dx,dy,x,y,Lx,Ly,i,j,w,Umax);    

        % Model horizontal flows
        fx(j,i) = fx(j,i) - kStiff*( uX(j,i) - uX_Tar );
        

        
        % Model minimal vertical flow through top and bottom
        if ( y(j) < 0.025 ) || ( y(j) > 0.475)
            fy(j,i) = fy(j,i) - 2.5*kStiff*( uY(j,i) - 2.5*Umax * (tanh(20*t)) ); 
            %fy(j,i) = fy(j,i) - kStiff*( uY(j,i) - 2.5*Umax * (tanh(20*t)) ); 
        else
            %fy(j,i) = fy(j,i);
            fy(j,i) = fy(j,i) - kStiff*( uY(j,i) - 0 ); 
        end
        
    end
%end
fx_exts = fx;
fy_exts = fy;

% MIGHT NOT NEED THESE!
%fx_exts = fx/ds^2;
%fy_exts = fy/ds^2;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes the Target Velocity Profile (MODEL DEPENDENT)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [uX_Tar] = please_Give_Target_Velocity(t,dx,dy,xGrid,yGrid,Lx,Ly,i,j,w,Umax)

% t:     current time in simulation
% dx:    x-Grid spacing
% dy:    y-Grid spacing
% xGrid: vector of xPts in Eulerian grid
% yGrid: vector of yPts in Eulerian grid
% Lx:    x-Length of Eulerian Grid
% Ly:    y-Length of Eulerian Grid
% i:     ith component in x-Grid
% j:     jth component in y-Grid
% w:     width of Channel
% Umax:  maximum velocity

y = yGrid(j);  % y-Value considered

% if y<0.125
%     %uY_Tar = 0;                        % No external forces in y-direction
%     uX_Tar = 0.5*Umax * (tanh(20*t));  % Only external forces in x-direction
% elseif y<0.25
%     %uY_Tar = 0;                        % No external forces in y-direction
%     uX_Tar = 1*Umax * (tanh(20*t));    % Only external forces in x-direction
% elseif y<0.375
%     %uY_Tar = 0;                       % No external forces in y-direction
%     uX_Tar = 1.5*Umax * (tanh(20*t));  % Only external forces in x-direction
% else
%     %uY_Tar = 0;                      % No external forces in y-direction
%     uX_Tar = 2*Umax * (tanh(20*t));   % Only external forces in x-direction
% end

sF=0.25; % speed scale factor
tF=20;   % time-scale factor for ramping up speed

if y<0.125
    %uY_Tar = 0;                        % No external forces in y-direction
    uX_Tar = sF * Umax * (tanh(tF*t));  % Only external forces in x-direction
elseif y<0.25
    %uY_Tar = 0;                         % No external forces in y-direction
    uX_Tar = -sF * Umax * (tanh(tF*t)); % Only external forces in x-direction
elseif y<0.375
    %uY_Tar = 0;                       % No external forces in y-direction
    uX_Tar = sF * Umax * (tanh(tF*t)); % Only external forces in x-direction
else
    %uY_Tar = 0;                         % No external forces in y-direction
    uX_Tar = -sF * Umax * (tanh(tF*t));  % Only external forces in x-direction
end



