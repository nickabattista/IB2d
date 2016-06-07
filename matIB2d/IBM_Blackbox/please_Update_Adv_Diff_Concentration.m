%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Setting up advection-diffusion solver
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C = please_Update_Adv_Diff_Concentration(C,dt,dx,dy,uX,uY,k)

% C:     concentration 
% dt:    time-step
% dx,dy: spatial steps in x and y, respectively
% uX:    x-Component of Velocity
% uY:    y-Component of Velocity
% k:     diffusion coefficient

% Compute Necessary Derivatives (Note: these calculations could be parallalized)
Cx = give_Necessary_Derivative(C,dx,uX,'x');
Cy = give_Necessary_Derivative(C,dy,uY,'y'); 
% Cx = D(C,dx,'x'); % Cy = D(C,dy,'y');
Cxx = DD(C,dx,'x');
Cyy = DD(C,dy,'y');
    
% Update Concentration 
% C = C + dt * ( k*(Cxx+Cyy) - uX'.*Cx - uY'.*Cy );

C = C + dt * ( k*(Cxx+Cyy) - uX.*Cx - uY.*Cy );


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

        %left side of grid
        if signs(i,1) <= 0 
            C_z(i,1) =  ( C(i,2) - C(i,1) ) / (dz);
        else
            C_z(i,1) =  ( C(i,1) - C(i,len) ) / (dz);
        end

        %right side of grid
        if signs(len,1) <= 0 
            C_z(i,len) =  ( C(i,1) - C(i,len) ) / (dz);
        else
            C_z(i,len) =  ( C(i,len) - C(i,len-1) ) / (dz);
        end

    end
    %Standard Upwind 
    for i=1:len
        for j=2:len-1
            if signs(i,j) <= 0
                C_z(i,j) = ( C(i,j+1) - C(i,j) ) / (dz);
            else
                C_z(i,j) = ( C(i,j) - C(i,j-1) ) / (dz);
            end
        end
    end

    % Ends x-Direction calculation %
    
elseif strcmp(string,'y') 

    %For periodicity on ends w/ UPWIND
    for i=1:len

        %bottom of grid
        if signs(1,i) <= 0 
            C_z(1,i) =  ( C(2,i) - C(1,i) ) / (dz);
        else
            C_z(1,i) =  ( C(1,i) - C(len,i) ) / (dz);
        end

        %top side of grid
        if signs(len,1) <= 0 
            C_z(len,i) =  ( C(1,i) - C(len,i) ) / (dz);
        else
            C_z(len,i) =  ( C(len,i) - C(len-1,i) ) / (dz);
        end

    end
    
    %Standard Upwind
    for i=2:len-1
        for j=1:len
            if signs(i,j) <= 0
                C_z(i,j) = ( C(i+1,j) - C(i,j) ) / (dz);
            else
                C_z(i,j) = ( C(i,j) - C(i-1,j) ) / (dz);
            end
        end
    end

    % Ends y-Direction calculation %
    
else
        
    fprintf('\n\n\n ERROR IN FUNCTION D FOR COMPUTING 1ST DERIVATIVE\n');
    fprintf('Need to specify which desired derivative, x or y.\n\n\n'); 
       
end
    
clear signs; clear len;

