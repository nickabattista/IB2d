function [v_INTERP,ind_Start,ind_End] = FitzHugh_Nagumo_1d(dt_IBM,T_final)
%function FitzHugh_Nagumo_1d(dt_IBM,T_final)

%
% This script solves the FitzHugh-Nagumo Equations in 1d, which are 
% a reduced order model of the more complicated Hodgkin-Huxley Equations. 
%
% Author:  Nick Battista
% Created: 09/11/2015
% University: UNC-CH
%
% Equations:
% dv/dt = D*Laplacian(v) + v*(v-a)*(v-1) - w - I(t)
% dw/dt = eps*(v-gamma*w)
%
% Inputs:
% dt_IBM: time-step from immersed boundary code
% T_final: final time for immersed boundary simulation
%
% Variables & Parameters:
% v(x,t): membrane potential
% w(x,t): blocking mechanism
% D:      diffusion rate of potential
% a:      threshold potential
% gamma:  resetting rate
% eps:    strength of blocking
% I(t):   initial condition for applied activation
%


%USER-INPUT FROM MODEL:
numPtsAlongTube = 78;    % # of Lag. Muscle Pts Along Tube 
ind_Start = 155;         % first index for muscle between tube in spring count
ind_End = 232;           % last index for muscle between tube in spring count
dt_IBM = dt_IBM*2500;    % scale time-info
T_final_Ca = T_final;    % Stays on time-scale of IBM
T_final = T_final * 2500;% scale time-info


% Save Movie? 
PLAY_MOVIE = 0; % (1 for YES, 0 for NO)

% Parameters in model %
D = 10.0;       % Diffusion coefficient
a = 0.3;        % Threshold potential (Note: a=0.3 is traveling wave value, a=0.335 is interesting)
gamma = 1.0;    % Resetting rate (Note: large values give 'funky thick' traveling wave, gamma = 1.0 is desired)
eps = 0.001;    % Blocking strength (Note: eps = 0.001 is desired)
I_mag = 0.05;   % Activation strength

% Discretization/Simulation Parameters %
factor = 4;     % 4x the resolution of the FHN model than IBM for tube
N = factor*numPtsAlongTube-1; % # of discretized points for FHN model (finer resolution than IBM mesh)
L = 500;        % Length of domain, [0,L]
dx = L/N;       % Spatial Step
x = 0:dx:L;     % Computational Domain

% Temporal  Parameters %
Np = 10;              % Set the number of pulses
pulse = T_final/Np;   % determines the length of time between pulses.
numFHN = 10;          % # that relates # of time-steps of FHN to IBM
dt = dt_IBM / numFHN; % Time-step for FitzHugh-Nagumo (fraction of IBM time-step)
NT = T_final / dt;    % Total # of time-steps to be taken for FHN
NT_IBM = T_final / dt_IBM; % Total # of time-steps for IBM
i1 = 0.25;%0.475;     % fraction of total length where current starts
i2 = 0.3;%0.525;      % fraction of total length where current ends
dp = pulse/50;        % Set the duration of the current pulse
pulse_time = 0;       % pulse time is used to store the time that the next pulse of current will happen
IIapp=zeros(1,N+1);   % this vector holds the values of the applied current along the length of the neuron
dptime = T_final/100; % This sets the length of time frames that are saved to make a movie.

% Initialization %
v = zeros(1,N+1);
w = v;
t=0;
ptime = 0;       
store = 1;       % counter for storing data to match time-steps from IBM
tVec = 0:dt:T_final;
Nsteps = length(tVec);
vStore = zeros(NT_IBM,N+1); vStore(store,:) = v;
wStore = zeros(NT_IBM,N+1); wStore(store,:) = w;
IIappStore=vStore;
store = store+1; % update storage counter


% Compute Calcium-Dynamics a-priori for activation wave!
% (returns Ca: free calcium ions, Caf: bound to filaments Ca-ions)
fprintf('     --> Solving Calcium Dynamics Model\n');
[Ca,Caf] = Calcium_Dynamics(Nsteps,T_final_Ca);
fprintf('     --> Finished calculating Calcium Dynamics\n');



%
% **** % **** BEGIN SIMULATION! **** % **** %
%
for i=2:Nsteps;
    
     % Update the time
    t = t+dt;                        
    
    % Give Laplacian
    DD_v_p = give_Me_Laplacian(v,dx);  
    
    % Gives activation wave (either from prescribed or Calcium-Dynamics model)
    %[IIapp,pulse_time] = Iapp(pulse_time,i1,i2,I_mag,N,pulse,dp,t,IIapp);
    IIapp = Iapp_from_Calcium_Dynamics(i,i1,i2,N,Caf);
    
    % Update potential and blocking mechanism, using Forward Euler
    vN = v + dt * ( D*DD_v_p - v.*(v-a).*(v-1) - w + IIapp );
    wN = w + dt * ( eps*( v - gamma*w ) );
    
    % Update time-steps
    v = vN;
    w = wN;
    
    % Store time-step values
    if mod(i-1,numFHN) == 0
        vStore(store,:) = v;
        wStore(store,:) = w;
        IIappStore(store,:) = IIapp;
        store = store + 1;
    end
    
    % PLAY MOVIE?
    if PLAY_MOVIE == 1
        %This is used to determine if the current time step will be a frame in the movie
        if t > ptime,
            figure(1)
            plot(x, v);
            axis([0 L -0.5 1.5]);
            xlabel('Distance (x)');
            ylabel('Electropotenital (v)');
            ptime = ptime+dptime;
            fprintf('Time(s): %d\n',t);
            pause(0.01);
        end
    end
    
end %END TIME-STEPPING LOOP

% Compute electro-potential at associated Lagrangian Pts.
v_INTERP = compute_IBM_Potential_At_IBM_Lag_Pts(factor,numPtsAlongTube,vStore);
%IIapp_INTERP = compute_IBM_Potential_At_IBM_Lag_Pts(factor,numPtsAlongTube,IIappStore);


% TEST SOLUTION
%  pause();
%  x=1:1:78;
%  for i=1:10:length(v_INTERP(:,1))
%     plot(x,v_INTERP(i,:),'-'); 
%     %plot(x,IIapp_INTERP(i,:),'-');
%     axis([0 79 -0.5 1.5]);
%     pause(0.01);
%  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes the electro-potential, v(x,t), at the correct points
% along the IBM's Lagrangian Structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [v_INTERP] = compute_IBM_Potential_At_IBM_Lag_Pts(factor,numPtsAlongTube,vStore)

% v: row: electro-potential 
% factor: resolution of FHN:IBM 
% vStore: membrane potential along tube

v_INTERP = zeros(length(vStore(:,1)),numPtsAlongTube);

ct = 1;
for j=1:length(vStore(1,:))
   if ( mod(j-1,factor) == 0 )
       v_INTERP(:,ct) = vStore(:,j);
       ct = ct+1;
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: the CALCIUM inspired injection function, Iapp = activation wave 
% for system, and returns the activation signal between i1 and i2 in
% geometry along tube
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function IIapp = Iapp_from_Calcium_Dynamics(ith_step,i1,i2,N,Caf)

    % Resets activation to zero
    IIapp = zeros(1,N+1);
    
    % Activates the proper region based on Calcium-Dynamics
    for j=(floor(i1*N):floor(i2*N))
        coeff = 0.5;                 % Coefficient to scale activation wave accordingly
        IIapp(j) = coeff*Caf(ith_step); % Activation of bound Calcium to filaments
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: the injection function, Iapp = activation wave for system, and
% returns both the activation as well as updated pulse_time
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [app,pulse_time] = Iapp(pulse_time,i1,i2,I_mag,N,pulse,dp,t,app)


    %Check to see if there should be a pulse
    if t > (pulse_time),
        
        % Sets pulsing region to current amplitude of I_mag x\in[i1*N,i2*N]
        for j=(floor(i1*N):floor(i2*N)),
            app(j) = I_mag;  
        end
        
        % Checks if the pulse is over & then resets pulse_time to the next pulse time.
        if t > (pulse_time+dp),
            pulse_time = pulse_time+pulse;
        end
        
    else
        
        % Resets to no activation
        app = zeros(1,N+1);
    
    end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives Laplacian of the membrane potential, note: assumes
% periodicity and uses the 2nd order central differencing operator.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DD_v = give_Me_Laplacian(v,dx)

Npts = length(v);
DD_v = zeros(1,Npts);

for i=1:Npts
   if i==1
       %DD_v(i) = ( v(i+1) - 2*v(i) + v(end) ) / dx^2;
       DD_v(i) = ( v(i+1) - 2*v(i) + 0 ) / dx^2;
   elseif i == Npts
       %DD_v(i) = ( v(1) - 2*v(i) + v(i-1) ) / dx^2;
       DD_v(i) = ( 0 - 2*v(i) + v(i-1) ) / dx^2;
   else
       DD_v(i) = ( v(i+1) - 2*v(i) +  v(i-1) ) /dx^2;
   end

end

