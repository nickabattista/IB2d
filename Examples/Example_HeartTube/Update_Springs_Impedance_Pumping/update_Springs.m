function springs = update_Springs(dt,current_time,xLag,springs)

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

freq = 10;  % Frequency for pumping
d = 1.0;    % Diameter of Heart Tube

for i=(N-2)+10:(N-2)+20            % Loops over desired springs!
    
    RL_Vec(i) = d - abs( 0.9*d*sin( 2*pi*freq*current_time ) );

end


springs(:,4) = RL_Vec; % Update the springs_info

