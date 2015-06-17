%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes a discrete approx. to a 1D Dirac-delta function over a
% specified matrix, x, and spatial step-size, dx. It will have support in
% [x-2dx, x+2dx]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function delta = compute_delta_kernel(x,dx)

% x:  Values in which the delta function will be evaulated
% dx: Spatial step-size of grid

% Computes Dirac-delta Approximation.
RMAT = abs(x)/dx;

%Initialize delta
delta = RMAT;

%Loops over to find delta approximation
[row,col] = size(x);
for i=1:row
    for j=1:col
        
        r = RMAT(i,j);
        
        if r<1
            delta(i,j) = ( (3 - 2*r + sqrt(1 + 4*r - 4*r.*r) ) / (8*dx) );
        elseif ( (r<2) && (r>=1) )
            delta(i,j) = ( (5 - 2*r - sqrt(-7 + 12*r - 4*r.*r) ) / (8*dx) );
        end
        
    end
end


