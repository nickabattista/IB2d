%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION finds the indices on the Eulerian grid where the 1D Dirac-delta
% kernel is possibly non-zero is x-dimension.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function indices = give_1D_NonZero_Delta_Indices(lagPts_j, N, dx, supp)

%lagPts_j: matrix of lagrangian pts for specific coordinate, j= x or y.
%N:        # spatial resolution of Eulerian grid in each dimension
%dx:       Spatial step-size on Eulerian (fluid) grid
%supp:     Size of support of the Dirac-delta kernel (should be even)


% Finds the index of the lower left Eulerian pt. to Lagrangian pt..
ind_Aux = floor(lagPts_j/dx + 1);

% Get all the different x indices that must be considered.
indices = [];
for i=1:supp
    indices = [indices ind_Aux]; 
end
%
for i=1:supp
    indices(:,i) = indices(:,i) + -supp/2+1+(i-1); 
end

% Translate indices between {1,2,..,N}
indices = mod(indices-1,N) + 1;
