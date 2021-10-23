%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Returns vector of muscle activation forces for a single Lag.Pt.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Fm = give_Muscle_Activation(v,LF,LFO,SK,a,b,Fmax,current_time,xPt,xLag)

% current_time: current time in simulation (s)
% xLag: vector of all x-Lagrangian Pts

% Fm = a_f * Fmax * F1(Lf) * F2(Vf) 
%    = a_f * Fmax * exp( -( (Q-1)/SK )^2 ) * (1/P0)*(b*P0-a*v)/(v+b); Q = LF/LFO
%

% a_f: coefficient that triggers contraction (traveling square or Gaussian wave?)

% Length Tension Model Parameters (F1)
% Fmax: maximum isometric force produced at the optimum length of the muscle fibers
% LF:   length of the muscle fibers
% LFO:  length at which the muscle fibers exert their maximum tension
% SK:   constant specific for each muscle where SK > 0.28.


% Hill Model Parameters (F2)
% P0:   maximum load w/ NO contraction
% a:    
% b:    
% v:    velocity of muscle expansion/contraction


% Length Tension Model Parameters %
Q = LF/LFO;
F1 = exp( - ( (Q-1)/SK )^2 );

% Hill Model %
P0 = Fmax;   %Same as Fmax
F2 = (1/P0)*(b*P0-a*v)/(v+b);

% Get Activation Coefficient %
af_Val = give_Traveling_Triggering_Coefficient(current_time,xLag,xPt);

% Actually Compute Muscle Force %
%Fm = size( xLag );
Fm = af_Val*Fmax*F1*F2;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Returns value of Activation Trigger at specific time
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function af_Val = give_Traveling_Triggering_Coefficient(current_time,xLag,xPt)

% current_time: current_time in simulation
% xLag:         x-Lagrangian pts (both bottom and top of tube)
% xPt: x-Pt of interest

t = current_time;                % current time
freq = 3;                        % frequency of traveling wave down tube
t = rem(t,1/freq);               % Gives remainder after "modular arithmetic" ("fmod" in C++)

af_Val = sin( 2*pi*freq*t );   % Gives activation coefficient, e.g., btwn [0,1] 

if t < 1/(2*freq)
   af_Val = 1;%sin( pi*freq*t );   % Gives activation coefficient, e.g., btwn [0,1] 
else
   af_Val = -2;%sin( pi*freq*t );  % Gives activation coefficient, e.g., btwn [0,1];  
end



