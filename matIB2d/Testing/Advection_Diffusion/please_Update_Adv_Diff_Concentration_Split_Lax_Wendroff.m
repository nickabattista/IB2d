%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Setting up advection-diffusion solver
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C = please_Update_Adv_Diff_Concentration_Split_Lax_Wendroff(C,dt,dx,dy,uX,uY,k)

% C:     concentration 
% dt:    time-step
% dx,dy: spatial steps in x and y, respectively
% uX:    x-Component of Velocity
% uY:    y-Component of Velocity
% k:     diffusion coefficient

% Compute Necessary Derivatives for x-Advection 
Cx = give_Necessary_Derivative(C,dx,uX,'x');
Cxx = DD(C,dx,'x');

% Time-step #1 (give auxillary)
C = C + dt * ( k*(Cxx) - uX.*Cx + 0.5*dt*(uX.^2).*Cxx );

% Compute Necessary Derivatives for y-Advection 
Cy = give_Necessary_Derivative(C,dy,uY,'y'); 
Cyy = DD(C,dy,'y');
 
% Time-step #2 (give next iteration)
C = C + dt * ( k*(Cyy) - uY.*Cy + 0.5*dt*(uY.^2).*Cyy );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Computes derivative based on sign of Velocity, u, using UPWIND
% approach
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C_z = give_Necessary_Derivative(C,dz,uZ,string)

C_z = zeros(size(C));
len = length(uZ(:,1));
signs = sign(uZ);

if strcmp(string,'x')

    %For periodicity on ends w/ UPWIND
    for i=1:len

        % left side of grid
        C_z(i,1) =  ( C(i,2) - C(i,len) ) / (2*dz);
        
        % right side of grid
        C_z(i,len)= ( C(i,1) - C(i,len-1) ) / (2*dz);

    end
    
    %Standard Central Differencing
    for i=1:len
        for j=2:len-1          
                C_z(i,j) = ( C(i,j+1) - C(i,j-1) ) / (2*dz);
        end
    end

    % Ends x-Direction calculation %
    
elseif strcmp(string,'y') 

    %For periodicity on ends w/ UPWIND
    for i=1:len

        
        %bottom of grid
        C_z(1,i) =  ( C(2,i) - C(len,i) ) / (2*dz);

        %top of grid
        C_z(len,i) =  ( C(1,i) - C(len-1,i) ) / (2*dz);
        
    end
    
    %Standard Central Differencing
    for i=2:len-1
        for j=1:len
                C_z(i,j) = ( C(i+1,j) - C(i-1,j) ) / (2*dz);
        end
    end

    % Ends y-Direction calculation %
    
else
        
    fprintf('\n\n\n ERROR IN FUNCTION D FOR COMPUTING 1ST DERIVATIVE\n');
    fprintf('Need to specify which desired derivative, x or y.\n\n\n'); 
       
end
    
clear signs; clear len;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Finds CENTERED finite difference approximation to 2ND
% DERIVATIVE in z direction, specified by input and 'string' 
% Note: It automatically accounts for periodicity of the domain.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u_zz = DD(u,dz,string)

% u:      velocity 
% dz:     spatial step in "z"-direction
% string: specifies which 2ND derivative to take (to enforce periodicity)

len = length(u(:,1));

if strcmp(string,'x')

    %For periodicity on ends
    u_zz(:,1) =  ( u(:,2) - 2*u(:,1)   + u(:,len) )   / (dz^2);
    u_zz(:,len)= ( u(:,1) - 2*u(:,len) + u(:,len-1) ) / (dz^2);

    %Standard Upwind Scheme (Centered Difference)
    for j=2:len-1
        u_zz(:,j) = ( u(:,j+1) - 2*u(:,j) + u(:,j-1) ) / (dz^2);
    end

elseif strcmp(string,'y')

    %For periodicity on ends
    u_zz(1,:) =  ( u(2,:) - 2*u(1,:)   + u(len,:) )   / (dz^2);
    u_zz(len,:)= ( u(1,:) - 2*u(len,:) + u(len-1,:) ) / (dz^2);

    %Standard Upwind Scheme (Centered Difference)
    for j=2:len-1
        u_zz(j,:) = ( u(j+1,:) - 2*u(j,:) + u(j-1,:) ) / (dz^2);
    end

else
    
    fprintf('\n\n\n ERROR IN FUNCTION DD FOR COMPUTING 2ND DERIVATIVE\n');
    fprintf('Need to specify which desired derivative, x or y.\n\n\n');  
    
end



