%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Nicholas A. Battista
% Institution: The College of New Jersey (TCNJ)
% Created: March 20, 2018
% Date of Last Revision: March 22, 2018
%
% Written for use in MAT 331: Numerical Analysis 
%
% PURPOSE:  -To compute the error ssociated with applying Trapezoid Rule to 
%            two different integral cases:
%                (1) a periodic integrand on a periodic integration domain
%                (2) non-periodic integrand
%
%           -To demonstrate exponential convergence of Trap. Rule in Case (1)
%            above
%
% Inputs:   
%
% Returns:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Trap_Rule_Error()

% Number of subintervals to test
NVec = [1e0:1e0:9e0 1e1:1e1:9e1 1e2:1e2:9e2 1e3:1e3:9e3 1e4:1e4:9e4 1e5:1e5:9e5 1e6:1e6:9e6];

%
% Storage Initialization
ErrVecPeriodic = zeros(length(NVec));
ErrVecNonPeriodic = ErrVecPeriodic;

%
% Compute Integral Approximation for Different #'s of Subintervals
%
for i=1:length(NVec)
   
   flag = 1; % For Periodic
   ErrVecPeriodic(i) = Trap_Rule( NVec(i), 0, 1, flag);
   
   flag = 0; % For Non-Periodic
   ErrVecNonPeriodic(i) = Trap_Rule( NVec(i), 0, 1, flag);
   
end


%
% Print Convergence Rate Information to Screen
slopeNonPeriodic = give_Me_Slope( ErrVecNonPeriodic,40,50 );
fprintf('\nConvergence Rate for Periodic Functions Appears to be Exponential, e.g.,\n');
fprintf('     we see a linear relationship between log(error) and # of subintervals, N\n');
fprintf('\nConvergence Rate for Non-Periodic Functions Appears to be: %d\n\n',slopeNonPeriodic);

%
% PLOT THE ABSOLUTEL ERROR -> GIVE CONVERGENCE PLOTS!
%
fs = 18;
ms = 30;
figure(1) % PERIODIC CASE
semilogy(NVec(1:45),ErrVecPeriodic(1:45),'b.','MarkerSize',ms); hold on;
%loglog(NVec(1:45),ErrVecPeriodic(1:45),'b.','MarkerSize',ms); hold on;
title('Periodic Case');
ylabel('Absolute Error');
xlabel('Number of Sub-Intervals');
set(gca,'FontSize',fs);

figure(2) % NON-PERIODIC CASE
loglog(NVec,ErrVecNonPeriodic,'b.','MarkerSize',ms); hold on;
title('Non-Periodic Case');
ylabel('Absolute Error');
xlabel('Number of Sub-Intervals');
set(gca,'FontSize',fs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes the absolute error using the trapezoid rule for a
%           particular number of subintervals, N.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function err = Trap_Rule(N,a,b,flagCase)

% Uniformly spaced integration nodes
x = a:(b-a)/N:b;

% Compute trapezoid rule for periodic function
int = 0;
for i=1:N
   int = int + (b-a)/(2*N)*( f(x(i),flagCase) + f(x(i+1),flagCase) ); 
end

%
% Return Error for Particular Case
%
if flagCase == 1
    
    % PERIODIC CASE
    intExactPeriodic = 0.132214293037990;
    err = abs( int - intExactPeriodic );

else

    % NON-PERIODIC CASE
    intExactNonPeriodic = 0.455122322888408;
    err = abs( int - intExactNonPeriodic );

end

% PLOT!
%x=a:0.001:b;
%for i=1:length(x)
%    plot( x(i) , f(x(i)), '*'); hold on;
%end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: returns integrand function value
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = f(x,flagCase)

%
% Return Function Value for Particular Case
%
if flagCase == 1

    % PERIODIC: f(x) = (cos(2*pi*x))^2 / ( 1 + e^(sin(2*pi*x)) )^2;
    val = ( cos(2*pi*x) )^2 / ( 1 + exp( sin(2*pi*x) ) )^2;

else
    
    % NON-PERIODIC: f(x) = (x^2+3)(cos(2*pi*x))^2 / ( 1 + e^(sin(2*pi*x)) )^2;
    val = (x^2+3)*(cos(2*pi*x))^2 / ( 1 + exp(sin(2*pi*x)) )^2;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: give me slope of convergence plot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function m = give_Me_Slope( errVec, val1, val2 )

m = ( log( errVec(val1) ) - log( errVec(val2) ) ) / (log(val1) - log(val2));