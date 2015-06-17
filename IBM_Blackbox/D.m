%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Finds CENTERED finite difference approximation to 1ST
% Derivative in specified direction by input, dz, and 'string'. 
% Note: It automatically accounts for periodicity of the domain.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u_z = D(u,dz,string)

% u:      velocity 
% dz:     spatial step in "z"-direction
% string: specifies which 1ST derivative to take (to enforce periodicity)


len = length(u(:,1));

if strcmp(string,'x')

    %For periodicity on ends
    u_z(:,1) = ( u(:,2) - u(:,len) ) / (2*dz);
    u_z(:,len)= ( u(:,1) - u(:,len-1) ) / (2*dz);

    %Standard Centered Difference 
    for j=2:len-1
        u_z(:,j) = ( u(:,j+1) - u(:,j-1) ) / (2*dz);
    end

elseif strcmp(string,'y')
    
    %For periodicity on ends
    u_z(1,:) = ( u(2,:) - u(len,:) ) / (2*dz);
    u_z(len,:)= ( u(1,:) - u(len-1,:) ) / (2*dz);

    %Standard Centered Difference 
    for j=2:len-1
        u_z(j,:) = ( u(j+1,:) - u(j-1,:) ) / (2*dz);
    end
    
else
    
    fprintf('\n\n\n ERROR IN FUNCTION D FOR COMPUTING 1ST DERIVATIVE\n');
    fprintf('Need to specify which desired derivative, x or y.\n\n\n'); 
       
end


