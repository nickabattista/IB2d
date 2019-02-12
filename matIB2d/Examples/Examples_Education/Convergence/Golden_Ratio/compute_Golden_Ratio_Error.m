%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Nicholas A. Battista
% Institution: The College of New Jersey (TCNJ)
% Created: January 11, 2018
% Date of Last Revision: January 12, 2018
%
% Written for use in MAT 331: Numerical Analysis 
%
% PURPOSE:  -Computes Golden Ratio Approximation for a Certain Number of Terms 
%            in Fibonacci Sequence (e.g., line 25)
%
%           -To compute the error associated with that approximation
%
% Inputs:   
%
% Returns:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function compute_Golden_Ratio_Error()

% How many Fibonacci Terms do we want?
N = 40; 

% Gives the first N Fibonacci Numbers to us
Fib_Vec = compute_Fibonacci_Terms(N);

% Computes the first N-1 approximate Golden Ratio Values
for i=1:N-1
    GR_Vec(i) = Fib_Vec(i+1) / Fib_Vec(i);
end

% True value of Golden Ratio (given)
GR_True = ( 1 + sqrt(5) ) / 2;

% Compute the Golden Ratio Error for all approximate values
for i=1:N-1
    
    % Calculate the absolute error
    ErrVec(i) = abs( GR_True - GR_Vec(i) );
    
end

% Plot the error vs. number of terms in the sequence
fs = 18; % FontSize
ms = 32; % MarkerSize
numTerms = 2:1:N; % Need at least 2 terms to approximate Golden Ratio
semilogy(numTerms,ErrVec,'.','MarkerSize',ms); hold on;
xlabel('# of terms in Fib. Sequence');
ylabel('Absolute Error');
set(gca,'FontSize',fs)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: FUNCTION: Computes first N terms in the Fibonacci Sequence
%
% Inputs:  
%           N: number of terms desired
%
% Returns:  
%           Fib_Vec: vector containing Fibonacci Terms
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Fib_Vec = compute_Fibonacci_Terms(N)

% N: # of terms desired in Fibonacci sequence

% Initialize first two spots in storage vector for first two Fibonacci numbers
Fib_Vec(1) = 1; % e.g., F_1 = 1
Fib_Vec(2) = 1; % e.g., F_2 = 1

for i=3:N
    
   % Compute next term in Fibonacci Sequence, F_{N} = F_{N-1} + F_{N-2}, use "i" as the indexing variable for storage
   Fib_Vec(i) = Fib_Vec(i-1) + Fib_Vec(i-2); 
   
end

%We now have a vector of Fibonacci Values for [F_1,F_2,...,F_N]
