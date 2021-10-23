%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes the angular displacements of a damped pendulum under
%           gravity
%
% PURPOSE: to illustrate that changing the pendulum bob's radius results in
%          significant dynamical differences in the pendulum system
%
% Author: N.A. Battista
% Date: 12/13/2019
% Institution: The College of New Jersey
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function damped_pendulum_diff_r()

%
% Physical Parameters
%
m = 100;    % Mass of Pendulum bob
g=9.80665;  % Gravitational Acceleration;
L=0.2;      % Length of Pendulum

% 
% Radii of Bobs (here - all different)
%
rVec = 100*[0.0025 0.005 0.01 0.015 0.02 0.025];

%
% Damping Parameter (here - all the same)
%
bVec = 50*ones(1,6);

% 
% Temporal Parameters
%
t=0;        % initial time
tFinal = 24; % final-time
dt = 1e-3;  % time-step
ct=1;       % counter
tVec(ct)=t; % time vector storage

%
% Moment of Inertia
% 
IVec = m*(0.5*rVec.^2+L^2); 

% 
% Initial Angular Displacement
%
ang = -(pi/2-pi/5)*ones(1,6);
phi = zeros( size( ang ) );

%
% Perform the time-stepping, e.g., solve the ODEs
%
while t<=tFinal
    
    % increment the storage counter
    ct = ct+1;

    % increment time
    t = t+dt;
   
    % save current time
    tVec(ct) = t;
    
    % solve 1st order ODE for angular displacement
    ang(ct,:) = ang(ct-1,:) + dt * ( phi(ct-1,:) );
    
    % solve 1st order ODE for angular velocity
    phi(ct,:) = phi(ct-1,:) + dt * ( -m*g*L ./IVec .* sin( ang(ct-1,:) ) - bVec.*phi(ct-1,:)./IVec );
    
end

%
% Defining Colors for Plotting Data
%
color2 = [0.635 0.078 0.184];  % dark red
color5 = [0.9 0.525 0.098];    % orange
color6 = [0.929 0.840 0.1250]; % dark yellow
color7 = [0.7 0.7 0];          % 'ugly' green
color13 = [0 0.5 0.5];             % green

%
% Note discrepancy is a bit off due to small angle approximation: 
% sin(theta) ~ theta gives below period of oscillation
%
T1 = 2*pi*sqrt( IVec(2) / (m*g*L) );
%
% Actually computed value of period is:
T1 = 1.93;

%
% Make Figure of Angular Displacement vs. Non-dimensional Time (# of oscillations of 'r' case)
%
lw=6;
fs=18;
plot(tVec/T1,ang(:,1),'k-','LineWidth',lw); hold on;
plot(tVec/T1,ang(:,2),'b-','LineWidth',lw,'Color',color2); hold on;
plot(tVec/T1,ang(:,3),'g-','LineWidth',lw,'Color',color5); hold on;
plot(tVec/T1,ang(:,4),'r-','LineWidth',lw,'Color',color6); hold on;
plot(tVec/T1,ang(:,5),'c-','LineWidth',lw,'Color',color7); hold on;
plot(tVec/T1,ang(:,6),'c-','LineWidth',lw,'Color',color13); hold on;
plot(tVec/T1,ang(:,5),'c-','LineWidth',lw,'Color',color7); hold on;
plot(tVec/T1,ang(:,4),'r-','LineWidth',lw,'Color',color6); hold on;
plot(tVec/T1,ang(:,3),'g-','LineWidth',lw,'Color',color5); hold on;
plot(tVec/T1,ang(:,2),'b-','LineWidth',lw,'Color',color2); hold on;
plot(tVec/T1,ang(:,1),'k-','LineWidth',lw); hold on;
xlabel('Non-Dimensional Time (# periods of r-case)');
ylabel('Angular Displacement (Radians)');
axis([0 12 -1 1.5]);
leg=legend('$\frac{1}{2}$r','r','2r','3r','4r','5r','Orientation','horizontal','interpreter','latex');
set(gca,'FontSize',fs);
set(leg,'FontSize',fs+2);
leg.NumColumns=3;

