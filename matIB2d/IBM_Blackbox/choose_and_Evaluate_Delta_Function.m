%------------------------------------------------------------------------------------------
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
%------------------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: evaluates a selected Regularized Delta Function below
%
%       Inputs:  r: 'radius' in x or y direction from Lag Pt.
%               dz: grid resolution 
%
%       Output: delta: delta function value when evaluated at r
%
%       Note: dividing by 'dz' here saves a bit of time, e.g., will not 
%             need to divide an entire sparse matrix by 'dz' later.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function delta = choose_and_Evaluate_Delta_Function(r,dz)

        %--------------------------------------------
        % OLD DELTA FUNCTION APPROXIMATION
        %--------------------------------------------
        %if r<=2
        %   delta = 0.25*(1+cos(pi*r/2))/dz;
        %else
        %   delta = 0;
        %end
        
        %----------------------------------------------------------
        % PESKIN 4-PT DISCRETE DELTA FUNCTION (see Peskin 2002)
        %----------------------------------------------------------
        if r<1
           delta = ( (3 - 2*r + sqrt(1 + 4*r - 4*r*r) ) / (8*dz) );
        elseif ( (r<2) && (r>=1) )
           delta = ( (5 - 2*r - sqrt(-7 + 12*r - 4*r*r) ) / (8*dz) );
        else
           delta = 0;
        end