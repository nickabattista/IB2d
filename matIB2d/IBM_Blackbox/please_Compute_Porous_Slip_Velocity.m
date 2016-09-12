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
% One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting lengths, etc
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
% 
% If you would like us %to add a specific muscle model, please let Nick (nick.battista@unc.edu) know.
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes Porous Slip Velocity based on Darcy's Law,
%
%           U_porous = -alpha <F_Lag,hat{n}> / | d vec{X}_Lag /ds |
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Por_Mat,nX,nY] = please_Compute_Porous_Slip_Velocity(ds,xLag,yLag,porous_info,F_Lag)

% xLag: vector of x-Pts associated w/ x-Lagrangian pts
% yLag: vector of y-Pts associated w/ y-Lagrangian pts
% porous_info:  col 1: lag-ids for porous media
%               col 2: x-Lag pts for lag-ids 
%               col 3: y-Lag pts for lag-ids
%               col 4: porosity coefficient 

% # of porous media pts.
Np = length( porous_info(:,1) );

% Initialize storage
%Por_X = zeros(length(xLag),1);
%Por_Y = Por_X;

% Compute Lagrangian Derivatives
[xL_s,yL_s] = give_Me_Lagrangian_Derivatives(ds,Np,porous_info);

% Compute Normal Vector (unit normals)
[nX,nY,sqrtNorm] = give_Me_Lagrangian_Normal_Vectors(xL_s,yL_s);


% Compute Porous Slip Velocity
Up_X = - ( porous_info(:,4) ) .* F_Lag( porous_info(:,1) ,1).*nX ./ sqrtNorm;
Up_Y = - ( porous_info(:,4) ) .* F_Lag( porous_info(:,1) ,2).*nY ./ sqrtNorm;
Por_Mat = [Up_X Up_Y];

% Store porous slip velocities in appropriate vector for adding to current Lagrangian Velocity Computation
%Por_X( porous_info(:,1) ) = Up_X;
%Por_Y( porous_info(:,1) ) = Up_Y;
%Por_Mat = [Por_X Por_Y];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes Lagrangian Derivatives
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xL_s,yL_s] = give_Me_Lagrangian_Derivatives(ds,Np,porous_info)

xL = porous_info(:,2);   % x-Values for Porous points
yL = porous_info(:,3);   % y-Values for Porous points
c = porous_info(:,5);    % choice of where it falls on stencil

xL_s = zeros(Np,1);
yL_s = zeros(Np,1);


for i=1:Np
   if c(i) == -2
       xL_s(i) = ( -25/12*xL(i) + 4*xL(i+1) - 3*xL(i+2) + 4/3*xL(i+3) - 1/4*xL(i+4) ) / ds;
       yL_s(i) = ( -25/12*yL(i) + 4*yL(i+1) - 3*yL(i+2) + 4/3*yL(i+3) - 1/4*yL(i+4) ) / ds;
       
   elseif c(i) == -1
       xL_s(i) = ( -0.25*xL(i-1) - 5/6*xL(i) + 1.5*xL(i+1) - 0.5*xL(i+2) + 1/12*xL(i+3) ) / ds;
       yL_s(i) = ( -0.25*yL(i-1) - 5/6*yL(i) + 1.5*yL(i+1) - 0.5*yL(i+2) + 1/12*yL(i+3) ) / ds;

   
   elseif c(i) == 0
       xL_s(i) = ( 1/12*xL(i-2) - 2/3*xL(i-1) + 2/3*xL(i+1) - 1/12*xL(i+2) ) / ds;
       yL_s(i) = ( 1/12*yL(i-2) - 2/3*yL(i-1) + 2/3*yL(i+1) - 1/12*yL(i+2) ) / ds;

       
   elseif c(i) == 1
       xL_s(i) = ( -1/12*xL(i-3) + 0.5*xL(i-2) - 1.5*xL(i-1) + 5/6*xL(i) + 0.25*xL(i+1) ) / ds;
       yL_s(i) = ( -1/12*yL(i-3) + 0.5*yL(i-2) - 1.5*yL(i-1) + 5/6*yL(i) + 0.25*yL(i+1) ) / ds;
            
   elseif c(i) == 2
       xL_s(i) = ( 0.25*xL(i-4) - 4/3*xL(i-3) + 3*xL(i-2) - 4*xL(i-1) + 25/12*xL(i) ) / ds;
       yL_s(i) = ( 0.25*yL(i-4) - 4/3*yL(i-3) + 3*yL(i-2) - 4*yL(i-1) + 25/12*yL(i) ) / ds;
       
   else
      fprintf('\n\n');
      error('Error: Bad file format inside your .porous file!\n'); 
   end
end

% IF CLOSED STRUCTURE!
% for i=1:Np
%     if i==1
%        xL_s(1) = ( xL(2) - xL(end) ) / (2*ds); 
%        yL_s(1) = ( yL(2) - yL(end) ) / (2*ds);
%     elseif i<Np
%        xL_s(i) = ( xL(i+1) - xL(i-1) ) / (2*ds); 
%        yL_s(i) = ( yL(i+1) - yL(i-1) ) / (2*ds);
%     else
%        xL_s(i) = ( xL(1) - xL(end-1) ) / (2*ds); 
%        yL_s(i) = ( yL(1) - yL(end-1) ) / (2*ds);
%     end
% end

% Matrices for calculating derivatives
%AL2 = [1 1 1 1 1; 0 1 2 3 4; 0 1 4 9 16; 0 1 8 27 64; 0 1 16 3^4 4^4];
%AL1 = [1 1 1 1 1; -1 0 1 2 3; 1 0 1 4 9; -1 0 1 8 27; 1 0 1 16 81];
%A0 = [1 1 1 1 1; -2 -1 0 1 2; 4 1 0 1 4; -8 -1 0 1 8; 16 1 0 1 16];
%AR1 = [1 1 1 1 1; -3 -2 -1 0 1; 9 4 1 0 1; -27 -8 -1 0 1; 81 16 1 0 1];
%AR2 = [1 1 1 1 1; -4 -3 -2 -1 0; 16 9 4 1 0; -64 -27 -8 -1 0; 256 81 16 1 0];
%rhs = [0 1 0 0 0]';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % FUNCTION: compute Porous Connections
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function connects_Por = give_Me_Porous_Connections(ds,Np,porous_info)
% 
% xL = porous_info(:,2);
% yL = porous_info(:,3);
% 
% start = 1; 
% last = Np;
% prev = last;
% for i=1:Np
%     
%     x0 = xL(i); y0 = yL(i);        % central node value
%     xm1= xL(prev); ym1 = yL(prev); % 'previous' values
%     xp1 =xL(i+1); yp1=yL(i+1);     % 'next' values
%     
%     distL = sqrt( (x0-xm1)^2 + (y0-ym1)^2 );
%     distR = sqrt( (x0-xp1)^2 + (y0-yp1)^2 );
%     
%     
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes Lagrangian UNIT Normal Vectors
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nX,nY,sqrtN] = give_Me_Lagrangian_Normal_Vectors(xL_s,yL_s)

sqrtN = sqrt( (xL_s).^2 + (yL_s).^2 );

nX = ( yL_s ) ./ sqrtN;
nY = ( -xL_s) ./ sqrtN;