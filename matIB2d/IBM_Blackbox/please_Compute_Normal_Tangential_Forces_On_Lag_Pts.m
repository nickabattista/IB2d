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
%   4. Mass Points
%   5. Porous Points
%	6. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-Tension)")
%   7. 3-Element Hill Muscle Model
%
% One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting-lengths, etc
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
% 
% If you would like us to add a specific muscle model, please let Nick (nick.battista@unc.edu) know.
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes normal and tangential forces on Lagrangian Boundary
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F_Tan_Mag,F_Normal_Mag] = please_Compute_Normal_Tangential_Forces_On_Lag_Pts(lagPts,F_Lag)


Npts = length(lagPts);
X = lagPts(:,1);
Y = lagPts(:,2);

% Compute Normal/Tangential Vectors
[nX,nY,sqrtN] = give_Me_Lagrangian_Normal_Vectors(Npts,X,Y);
[tX,tY,~] = give_Me_Lagrangian_Tangent_Vectors(Npts,nX,nY,sqrtN);

% Project Force Data onto Normal / Tangent Vectors
[F_Tan,F_Normal] = give_Tangent_and_Normal_Force_Projections(Npts,F_Lag(:,1),F_Lag(:,2),nX,nY,tX,tY);
    
% Compute Colormap Force Magnitude Scalings
[F_Tan_Mag,F_Normal_Mag] = give_Force_Magnitude_Scalings(Npts,F_Tan,F_Normal);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes Lagrangian Derivatives for Normal/Tangential Vector Computation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xL_s,yL_s] = give_Me_Lagrangian_Derivatives(ds,Npts,X,Y)

xL = X;
yL = Y;

xL_s = zeros(Npts,1);
yL_s = zeros(Npts,1);

for i=1:Npts
    if i==1
       xL_s(1) = ( xL(2) - xL(end) ) / (2*ds); 
       yL_s(1) = ( yL(2) - yL(end) ) / (2*ds);
    elseif i<Npts
       xL_s(i) = ( xL(i+1) - xL(i-1) ) / (2*ds); 
       yL_s(i) = ( yL(i+1) - yL(i-1) ) / (2*ds);
    else
       xL_s(i) = ( xL(1) - xL(end-1) ) / (2*ds); 
       yL_s(i) = ( yL(1) - yL(end-1) ) / (2*ds);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes Lagrangian UNIT Normal Vectors
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nX,nY,sqrtN] = give_Me_Lagrangian_Normal_Vectors(Npts,X,Y)

% Compute Lagrangian Spacing
ds = sqrt( ( X(3)-X(4) )^2 + ( Y(3)-Y(4) )^2 );

% Gives Lagrangian Derivatives
[xL_s,yL_s] = give_Me_Lagrangian_Derivatives(ds,Npts,X,Y);

sqrtN = sqrt( (xL_s).^2 + (yL_s).^2 );

nX = ( yL_s ) ./ sqrtN;
nY = ( -xL_s) ./ sqrtN;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes Lagrangian UNIT Tangent Vectors
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tX,tY,sqrtN] = give_Me_Lagrangian_Tangent_Vectors(Npts,nX,nY,sqrtN)

% Allocate storage
tX = zeros( size(nX) );
tY = zeros( size(nY) );

% Rotate normal vectors to get tangent vectors
ang = -pi/2; % Rotate CW by 90 degrees
for i=1:Npts
    tX(i) = nX(i)*cos(ang) - nY(i)*sin(ang);
    tY(i) = nX(i)*sin(ang) + nY(i)*cos(ang);
end

% Testing
%for i=1:Npts
%   test(i) = ( nX(i)*tX(i) + nY(i)*tY(i) ); 
%   test2(i) = sqrt( tX(i)*tX(i) + tY(i)*tY(i) );
%end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes force vector projections onto the tangent and normal
%           vectors! 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [F_Tan,F_Normal] = give_Tangent_and_Normal_Force_Projections(Npts,Fx,Fy,nX,nY,tX,tY)

% Allocate Storage
F_Tan = zeros(Npts,2); F_Normal = F_Tan;

for i=1:Npts

    % Compute dot product between force vector and tangent vector
    tanVec_dotProd = ( Fx(i)*tX(i) + Fy(i)*tY(i) ) / sqrt( tX(i)*tX(i) + tY(i)*tY(i) );
    F_Tan(i,1) = tanVec_dotProd * ( tX(i) );
    F_Tan(i,2) = tanVec_dotProd * ( tY(i) );
    
    % Compute dot product between force vector and normal vector
    normalVec_dotProd = ( Fx(i)*nX(i) + Fy(i)*nY(i) ) / sqrt( nX(i)*nX(i) + nY(i)*nY(i) );
    F_Normal(i,1) = normalVec_dotProd * ( nX(i) );
    F_Normal(i,2) = normalVec_dotProd * ( nY(i) );
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: scales the force matrices by desired percentiles of each in magnitude 
%           for colormap scalings 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MagTan,MagNormal] = give_Force_Magnitude_Scalings(Npts,F_Tan,F_Normal)

% Allocate Storage
MagTan = zeros( length( F_Tan(:,1) ),1 );
MagNormal = MagTan;

% Find the magnitude of the force in each 
for i=1:Npts
    MagTan(i,1) = sqrt( F_Tan(i,1)^2  + F_Tan(i,2)^2 );
    MagNormal(i,1)= sqrt( F_Normal(i,1)^2 + F_Normal(i,2)^2);
end

% % Finds Percentiles for Forces for Scaling
% prc90_T = prctile(MagTan,90);
% prc90_N = prctile(MagNormal,90);
% prc10_T = prctile(MagTan,10);
% prc10_N = prctile(MagNormal,10);
% 
% % "Scales the forces" via if-elseif statements by desired percentiles. 
% for i=1:Npts
%     
%     mT = MagTan(i);
%     mN = MagNormal(i);
%     
%     if mT >= prc90_T
%         MagTan(i) = prc10_T;
%     elseif mT <= prc10_T
%         MagTan(i) = prc10_T;
%     end
%     
%     if mN >= prc90_N
%         MagNormal(i) = prc10_N;
%     elseif mN <= prc10_N
%         MagNormal(i) = prc10_N;
%     end
%     
% end
