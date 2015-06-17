%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Returns vector of muscle activation forces for all xLag Pts.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Fm = give_Muscle_Activation(current_time,xLag,yLag)

% current_time: current time in simulation (s)
% xLag: vector of all x-Lagrangian Pts

% Fm = a_f * Fmax * F1(Lf) * F2(Vf) 
%    = a_f * F_IO*exp( -( (Q-1)/SK )^2 ) * (1/P0)*(b*P0-a*v)/(v+b); Q = LF/LFO
%

% a_f: coefficient that triggers contraction (traveling square or Gaussian wave?)

% Length Tension Model Parameters (F1)
% F_IO: maximum isometric force produced at the optimum length of the muscle fibers
% LF:   length of the muscle fibers
% LFO:  length at which the muscle fibers exert their maximum tension
% SK:   constant specific for each muscle where SK > 0.28.

% Hill Model Parameters (F2)
% P0:   maximum load w/ NO contraction
% a:    
% b:    
% v:

D = 1.0;     % diameter (width) of HeartTube
Fmax = 5e2;  % Max Force

% Length Tension Model Parameters %
LFO = D;
LF = D; %*********should be the current length of the muscle****************
Q = LF/LFO;
SK = 0.3;
F1 = exp( - ( (Q-1)/SK )^2 );

% Hill Model %
v = 0; %**************should be the current velocity of the muscle*********
a = 0.25;
b = 4.0;

% OTHER PIECES TO COMBINED MODEL %
%Fmax = 1.0; %TRY MULTIPLE VALUES!
P0 = Fmax;   %Same as Fmax
F2 = (1/P0)*(b*P0-a*v)/(v+b);

af_Vec = give_Traveling_Triggering_Coefficient(current_time,xLag);
Fm = zeros(length(xLag),1);      %initialize storage
for i=1:length(xLag)
   Fm(i,1) = af_Vec(i)*Fmax*F1*F2; 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Returns value of Activation Trigger at specific time
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function af_Vec = give_Traveling_Triggering_Coefficient(current_time,xLag)

% current_time: current_time in simulation
% xLag:         x-Lagrangian pts (both bottom and top of tube)

t = current_time;               % current time
freq = 1.0;                     % frequency of traveling wave down tube
Ltube = xLag(end/2) - xLag(1);  % length of heart tube
buff = 0.5;                     % percent-buff/2 on each end of heart tube
L_AR = Ltube*(1-buff);          % length of the actual activation region of tube
SqWidth = L_AR/12;              % width of the square traveling activation wave 
v = (L_AR-SqWidth) * freq;      % traveling wave velocity

t = rem(t,freq);                 % gives remainder after division ("fmod" in C++)

xL = xLag(1)+ (buff/2)*L_AR;  %+ (v*t);
xR = xLag(1)+ (buff/2)*L_AR  +L_AR; % SqWidth + (v*t);


af_Vec = zeros(1,length(xLag)); %initialize storage for activation coefficients
for i=1:length(xLag)
    x = xLag(i);

    if ( ( x >= xL ) && ( x <= xR ) )
        af_Vec(i) = (1/2)*(sin(2*pi*freq*current_time-pi/2)+1);
    else
        af_Vec(i) = 0.0;
    end
end




