%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Updates the concentration via a WENO Advection/Diffusion Scheme
%
%   Author: Matea Santiago
%   Date: February 2021
%   Institution (when created): UC Merced
%
%       Returns: Concentration on Eulerian Grid
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[C]=Adv_Diff_source(C,U,V,XEx,YEx,dt,D,LC2c1,LC2c2,permLC2,Nx,Ny,dx,dy,f)

% Half step of third order WENO scheme to solve advection equation
[C]=WENO_3O(C,C,U,V,XEx,YEx,Nx,Ny,dx,dy,dt/2);

% Full step of Crank Nicolson to solve Diffusion Equation
[C]=CN(C,D,LC2c1,LC2c2,permLC2,dx,dy,dt,Nx,Ny,f);

% Half step of third order WENO scheme to solve advection equation
[C]=WENO_3O(C,C,U,V,XEx,YEx,Nx,Ny,dx,dy,dt/2);

end
