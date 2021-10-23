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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: updates Spring parameters (either resting length of spring stiffness)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function springs = update_Springs(dt,current_time,xLag,yLag,springs)

% dt: time-step
% current_time: current time of simulation
% xLag: Vector of all xLagrangian Pts        
% springs: col 1: starting spring pt (by lag. discretization)
%          col 2: ending spring pt. (by lag. discretization)
%          col 3: spring stiffness
%          col 4: spring resting lengths

s_1 = springs(:,1);    % MASTER Node Indices
s_2 = springs(:,2);    % SLAVE Node Indices
k_Vec = springs(:,3);  % Springs Stiffnesses Vector
RL_Vec= springs(:,4);  % Resting Lengths Vectors

N = length(xLag);                   % Gives total number of Lagrangian pts!
N_springs = length( springs(:,1) ); % Gives total number of springs!

d = 1.0;  
L  = 2.0;
freq = 4;
v = L*freq;
t = rem(current_time, 1/freq);

%xL = xLag(20);% + v*current_time;
%xR = xLag(60);% + v*current_time;

% xL = xLag(20) + v*t;
% xR = xLag(25) + v*t;
% for i=1:N/2
%    xPt = xLag(i);
%    if ( ( xPt >= xL ) && ( xPt <= xR ) )
%         ii = (N-2) + i;
%         RL_Vec(ii) = d - abs( 0.9*d*sin( 2*pi*freq/2*current_time ) );
%    else
%         ii = (N-2) + i;
%         RL_Vec(ii) = d;
%    end
% end

xM = xLag(N/4);
for i=(N-2)+2:N_springs-1;
    
   id = i - (N-2); 
   x = xLag(id);
   delta = sin( (2*pi/L)*( (x - xM) - v*t ) );
   
   if delta >= 0
      RL_Vec(i) = d - 0.95*delta; 
   else
      RL_Vec(i) = d; 
   end
    
end

%pause();
springs(:,4) = RL_Vec; % Update the springs_info

