function [Ca,Caf] = Calcium_Dynamics(Nsteps,t_Final)
%function Calcium_Dynamics(Nsteps,t_Final)

%
% This script uses the  law of mass action to model the rates that calcium
% ions are bound and then released from the actin filaments and the 
% sarcoplasmic reticulum (SR). 
%
% This model, the Williams et al. 2008 model, assumes that the stimulus from 
% the motor neurons is either on or off, and the dynamics of the action 
% potentials in the motor neuron,  the release of acetylcholine in the synapse, 
% the binding of acetylcholine to receptors on the muscle cell membrane, and the 
% resulting depolarization of the muscle cell are neglected. 
%
%
% Author: Nicholas A. Battista
% Written: February 14, 2016
% University: UNC-CH
%
%
% Model Eqns:
% d(Ca)/dt = (k4*Caf-k3*Ca)*(1-Caf)+k1*(C-Ca-Caf)        [STIMULUS ON!]
% or
% d(Ca)/dt = (k4*Caf-k3*Ca)*(1-Caf)+k2*(Ca*(C-S-Ca-Caf)) [STIMULUS OFF!]
% and
% d(Caf)/dt = -(k4*Caf-k3*Ca)*(1-Caf)
%
%
% Variables & Parameters:
% (BEFORE SCALING)
% c:  free calcium ions
% s:  unbound SR calcium-binding sites
% cs: calcium-bound SR sites
% f:  unbound contractile filament calcium-binding sites
% cf: calcium-bound filament sites
% k1-k4: rate constants of binding/release
% 
% (CONSERVED QUANTITIES-> # of calcium ions, SR binding sites, and filament 
% binding sites per liter are conserved)
% cs+c+cf = CT (# of calcium ions)
% cs+s = ST    (# of SR binding sites)
% cf+f = FT    (# of filament sites)
%
% (AFTER SCALING)
% Caf = cf/FT (Caf <= 1, equality when all filaments are bound)
% Ca = c/FT   (Ca  <= C)
% C = CT/FT   (C is large enough for filament binding sites to be saturated
%              during tetanic stimulation)
% S = ST/FT   (S is large enough to reduce free calcium to a negligible
%              amount during rest)
%
% Units:
% length: mm = millimeters  (milli=10^-3)
% time:   s  = seconds
% force:  mN = milliNewtons (milli=10^-3; Newton=Kg m/s^2)
% velocity:  mm/s = millimeters per second
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Binding / Release Parameters
k1 = 9;                % Rate constant, Ca2+ binding in SR, s^-1
k2 = 50;               % Rate constant, Ca2+ release from SR, s^-1
k3 = 40;               % Rate constant, Ca2+ binding to filaments, s^-1
k4 = 19.4;             % Rate constant, Ca2+ release from filaments, s^-1

% Dimensionless Concentrations / Conserved Quantities
C = 200;               % total dimensionless Ca2+ concentration
S = 200;               % total dimensionless concentrations of sarcoplasmic-reticular binding sites


% Temporal Parameters for Ca-Dynamics / IBM Simulation
dt = t_Final/Nsteps;   %compute Ca-dynamics time-step

% Temporal Parameters for Activation
period = 0.5;          % period of Ca-wave
act_period = 0.005;    % actuation length of Ca-wave

% Initialization / Initial Conditions
Ca = zeros(Nsteps,1);
Caf = zeros(Nsteps,1);
Ca(1) = 0;             % inital free-calcium
Caf(1) = 0;            % initial bound calcium to filaments

%
% BEGIN TIME-STEPPING ROUTINE!
%

for i=2:Nsteps;
    
    % Calculate all RHS of DEQs at i-th time step
    uptake_Ca = k2*Ca(i)*(C-S-Ca(i)-Caf(i));       % uptake of calcium always happens
    
    % if stimulus is on -> free Ca is released by SR!
    if rem(dt*i,period)<act_period
        uptake_Ca = uptake_Ca + k1*(C-Ca(i)-Caf(i));
    end

    % Compute RHS of ODEs (removed weird non-biological term)
     RHS_Ca = k4*Caf(i) - k3*Ca(i)*(1-Caf(i)) + uptake_Ca;
     RHS_Caf = -k4*Caf(i) + k3*Ca(i)*(1-Caf(i));
    
    % Time-Step w/ Euler Method
    Ca(i+1) = Ca(i) + RHS_Ca*dt;
    Caf(i+1) = Caf(i) + RHS_Caf*dt;
    
end

figure(1)
time = 0:dt:t_Final;           %vector of times to solve ODEs
plot(time,Ca,time,Caf,'linewidth',2)
xlabel('Dimensionless Time');
ylabel('Dimensionless Concentration');
leg=legend(' [Ca]: free Ca ions', ' [Caf]: Ca bound filament sites');
maxY = 1.075*max(max(Ca,Caf));
axis([0 1.7 0 maxY] );
set(gca,'fontsize',18)
set(leg,'fontsize',18)
pause();
