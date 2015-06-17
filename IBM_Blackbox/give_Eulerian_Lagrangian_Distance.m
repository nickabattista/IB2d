%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION distance between Eulerian grid data, x, and Lagrangian grid data, y, 
%          at specifed pts typically and makes sure the distance are [0,L] accordingly.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function distance = give_Eulerian_Lagrangian_Distance(x, y, L)

%x,y: two matrices that you find the distance between (x-typically Eulerian data, y-typically Lagrangian data)
%L:   length of domain, i.e., [0,L]

[row,col] = size(x);
distance = abs( x - y );
for i=1:row
    for j=1:col
        distance(i,j) = min( distance(i,j), L-distance(i,j) ); %Note: need to make sure that taking correct value
    end
end

