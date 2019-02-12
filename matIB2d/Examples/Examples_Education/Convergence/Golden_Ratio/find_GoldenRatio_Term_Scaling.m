%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Nicholas A. Battista
% Institution: The College of New Jersey (TCNJ)
% Created: January 11, 2018
% Date of Last Revision: January 12, 2018
%
% Written for use in MAT 331: Numerical Analysis 
%
% PURPOSE:  -Finds scaling relation for # of terms in Fibonacci Sequence to
%           approximate the Golden Ratio to a certain error tolerance
%
%           -To demonstrate quick convergence to Golden Ratio to 1e-15
%
%           -To illustrate machine precision limitations numerically
%
% Inputs:   
%
% Returns:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function find_GoldenRatio_Term_Scaling()
                              

% Decide what error tolerance you want to test
errVec = [1e-1 1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10 1e-11 1e-12 1e-13 1e-14 1e-15 1e-16 1e-17 1e-18 1e-19 1e-20]; 


% Loop over all possible error tolerance and compute # of terms necessary for specified error
for i=1:length(errVec)
    
    % define the x as index i from the vector of possible values
    err = errVec(i);      
    
    % use function from earlier to find number of necessary terms in Taylor Series and save it to a new vector
    num_vector(i) = Golden_Ratio_Error(err);
    
end

%
% plots # of Fib.Seq. Necessary for a particular Error Tolerance with SemiLog plot on horizontal axis
%
figure(1)
fs = 18; % FontSize
ms = 32; % MarkerSize
semilogx(errVec,num_vector,'b.','MarkerSize',ms); hold on;
xlabel('Error Tolerances');
ylabel('# of Necessary Terms in Fib. Seq.');
set(gca,'FontSize',fs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: for a provided error tolerance, computes # of terms necessary
% in Fibonacci Sequence for approximate Golden Ratio value.
%
% Inputs:       
%          err_tol: error tolerance specified    
%
% Returns:     
%          num: # of terms necessary in Fibonacci Sequence
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function num = Golden_Ratio_Error(err_tol)


% Initialize Fibonacci Storage vector
Fib_Vec(1) = 1;
Fib_Vec(2) = 1;

% Initialization For While Loop
error = 1;   % could be any value > err_tol
i = 2;       % initialize counter for while loop
%
while error > err_tol
    
    % Iterate counter (so first round is F_3)
    i=i+1;
    
    % Compute next term in Fibonacci sequence
    Fib_Vec(i) = Fib_Vec(i-1)+Fib_Vec(i-2);
    
    % Compute Approximate "Golden Ratio"
    GR = Fib_Vec(i)/Fib_Vec(i-1);
    
    % Compute The Error
    error = abs( GR - ( 1+sqrt(5) )/2 );
    
end

% Save value of 
num = i;