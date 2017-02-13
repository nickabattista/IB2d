%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Setting up advection-diffusion solver
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C] = please_Update_Adv_Diff_Concentration_Flux_Limiter(C,dt,dx,dy,uX,uY,k)

% C:     concentration 
% dt:    time-step
% dx,dy: spatial steps in x and y, respectively
% uX:    x-Component of Velocity
% uY:    y-Component of Velocity
% k:     diffusion coefficient

% Compute Fluxes (Note: these calculations could be parallalized)
selection = 'superbee'; % choices: 'superbee' or 'vanLeer'

% compute upwind/downwind type approximations
cForwardx = [C(:,2:end) C(:,1)] - C(:,1:end);      % Gives Forward Euler
cBackwardx = C(:,1:end) - [C(:,end) C(:,1:end-1)]; % Gives Backward Euler

% compute fluxes in x-direction
F = give_Necessary_Fluxes(C,dx,uX,'x',selection,dt,cForwardx,cBackwardx); % Fluxes in x

% Compute 2nd Derivative Term in x
Cxx = DD(C,dx,'x');

% auxillary concentration
C = C + F + dt*( k*Cxx );
clear cForwardx cBackwardx Cxx F;

% compute upwind/downwind type approximations in y
cForwardy = [C(2:end,:); C(1,:)] - C(1:end,:);      % Gives Forward Euler
cBackwardy = C(1:end,:) - [C(end,:); C(1:end-1,:)]; % Gives Backward Euler

% compute fluxes in y-direction
F = give_Necessary_Fluxes(C,dy,uY,'y',selection,dt,cForwardy,cBackwardy); % Fluxes in y

% Compute 2nd Derivative Term in y
Cyy = DD(C,dy,'y');

% compute new concentration iterate
C = C + F + dt*( k*Cyy );
clear cForwardy cBackwardy Cyy F;






    
% UPWIND
%C = C - dt * ( diffX + diffY ) + dt*( k*laplacian_C );



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Computes derivative based on sign of Velocity, u, using UPWIND
% approach
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = give_Necessary_Fluxes(C,dz,uZ,string,selection,dt,cFz,cBz)

F = zeros(size(C));
len = length(uZ(:,1));
signs = sign(uZ);

if strcmp(string,'x')

    %For periodicity on ends w/ UPWIND
    for i=1:len

        %left side of grid
        if signs(i,1) <= 0         % u < 0
            lamX = uZ(i,1)*dt/dz;
            thetaP =   ( C(i,3) - C(i,2) ) / ( C(i,2) - C(i,1) );
            thetaM =   ( C(i,2) - C(i,1) ) / ( C(i,1) - C(i,len) );
            phi_P = please_Give_Flux_Limiter(thetaP,selection);
            phi_M = please_Give_Flux_Limiter(thetaM,selection);

            F(i,1) =  -lamX*cFz(i,1) - 0.5*lamX*(lamX-1)*( phi_M*cBz(i,1) - phi_P*cFz(i,1) );
        
        else
            
            lamX = uZ(i,1)*dt/dz;   % u > 0
            thetaP =   ( C(i,1) - C(i,len) ) / ( C(i,2) - C(i,1) );
            thetaM =   ( C(i,len) - C(i,len-1) ) / ( C(i,1) - C(i,len) );
            phi_P = please_Give_Flux_Limiter(thetaP,selection);
            phi_M = please_Give_Flux_Limiter(thetaM,selection);
            
            F(i,1) =  -lamX*cBz(i,1) - 0.5*lamX*(lamX-1)*( phi_P*cFz(i,1) - phi_M*cBz(i,1) );
            
        end

        %right side of grid
        if signs(i,len) <= 0    % u < 0
            
            lamX = uZ(i,len)*dt/dz;
            thetaP = ( C(i,2) - C(i,1) ) / ( C(i,1) - C(i,len) );
            thetaM = ( C(i,1) - C(i,len) ) / ( C(i,len) - C(i,len-1) );
            phi_P = please_Give_Flux_Limiter(thetaP,selection);
            phi_M = please_Give_Flux_Limiter(thetaM,selection);

            F(i,len) =  -lamX*cFz(i,len) - 0.5*lamX*(lamX-1)*( phi_M*cBz(i,len) - phi_P*cFz(i,len) );
            
        else                   % u > 0
            
            lamX = uZ(i,len)*dt/dz;
            thetaP = ( C(i,len) - C(i,len-1) ) / ( C(i,1) - C(i,len) );
            thetaM = ( C(i,len-1) - C(i,len-2) ) / ( C(i,len) - C(i,len-1) );
            phi_P = please_Give_Flux_Limiter(thetaP,selection);
            phi_M = please_Give_Flux_Limiter(thetaM,selection);
            
            F(i,len) =  -lamX*cBz(i,len) - 0.5*lamX*(lamX-1)*( phi_P*cFz(i,len) - phi_M*cBz(i,len) );
        end

    end
    
    %Standard Godunov Flux Limiter 
    for i=1:len
        for j=2:len-1
            if signs(i,j) <= 0 % u < 0
                
                if j==len-1
                    thetaP = ( C(i,1) - C(i,j+1) ) / ( C(i,j+1) - C(i,j) );
                    thetaM = ( C(i,j+1) - C(i,j) ) / ( C(i,j) - C(i,j-1) );
                else
                    thetaP = ( C(i,j+2) - C(i,j+1) ) / ( C(i,j+1) - C(i,j) );
                    thetaM = ( C(i,j+1) - C(i,j) ) / ( C(i,j) - C(i,j-1) );
                end
                
                lamX = uZ(i,j)*dt/dz;

                phi_P = please_Give_Flux_Limiter(thetaP,selection);
                phi_M = please_Give_Flux_Limiter(thetaM,selection);

                F(i,j) = -lamX*cFz(i,j) - 0.5*lamX*(lamX-1)*( phi_M*cBz(i,j) - phi_P*cFz(i,j) );
                
            else               % u > 0
                
                if j==2
                    thetaP = ( C(i,j) - C(i,j-1) ) / ( C(i,j+1) - C(i,j) );
                    thetaM = ( C(i,j-1) - C(i,len) ) / ( C(i,j) - C(i,j-1) );
                else
                    thetaP = ( C(i,j) - C(i,j-1) ) / ( C(i,j+1) - C(i,j) );
                    thetaM = ( C(i,j-1) - C(i,j-2) ) / ( C(i,j) - C(i,j-1) );
                end
                
                lamX = uZ(i,j)*dt/dz;

                phi_P = please_Give_Flux_Limiter(thetaP,selection);
                phi_M = please_Give_Flux_Limiter(thetaM,selection);
                
                F(i,j) = -lamX*cFz(i,j) - 0.5*lamX*(lamX-1)*( phi_P*cFz(i,j) - phi_M*cBz(i,j) );
            end
        end
    end

    % Ends x-Direction calculation %
    
elseif strcmp(string,'y') 

    %For periodicity on ends w/ UPWIND
    for i=1:len

        %bottom of grid
        if signs(1,i) <= 0        % u < 0
            lamY = uZ(1,i)*dt/dz;
            thetaP =   ( C(3,i) - C(2,i) ) / ( C(2,i) - C(1,i) );
            thetaM =   ( C(2,i) - C(1,i) ) / ( C(1,i) - C(len,i) );
            phi_P = please_Give_Flux_Limiter(thetaP,selection);
            phi_M = please_Give_Flux_Limiter(thetaM,selection);

            F(1,i) =  -lamY*cFz(1,i) - 0.5*lamY*(lamY-1)*( phi_M*cBz(1,i) - phi_P*cFz(1,i) );
        
        else                      % u > 0
            
            lamY = uZ(1,i)*dt/dz;
            thetaP =   ( C(1,i) - C(len,i) ) / ( C(2,i) - C(1,i) );
            thetaM =   ( C(len,i) - C(len-1,i) ) / ( C(1,i) - C(len,i) );
            phi_P = please_Give_Flux_Limiter(thetaP,selection);
            phi_M = please_Give_Flux_Limiter(thetaM,selection);
            
            F(1,i) =  -lamY*cBz(1,i) - 0.5*lamY*(lamY-1)*( phi_P*cFz(1,i) - phi_M*cBz(1,i) );
            
        end
        

        %top  of grid
        if signs(len,i) <= 0     % u < 0
            
            lamY = uZ(len,i)*dt/dz;
            thetaP = ( C(2,i) - C(1,i) ) / ( C(1,i) - C(len,i) );
            thetaM = ( C(1,i) - C(len,i) ) / ( C(len,i) - C(len-1,i) );
            phi_P = please_Give_Flux_Limiter(thetaP,selection);
            phi_M = please_Give_Flux_Limiter(thetaM,selection);

            F(len,i) =  -lamY*cFz(len,i) - 0.5*lamY*(lamY-1)*( phi_M*cBz(len,i) - phi_P*cFz(len,i) );
            
        else                    % u > 0
            
            lamY = uZ(len,i)*dt/dz;
            thetaP = ( C(len,i) - C(len-1,i) ) / ( C(1,i) - C(len,i) );
            thetaM = ( C(len-1,i) - C(len-2,i) ) / ( C(len,i) - C(len-1,i) );
            phi_P = please_Give_Flux_Limiter(thetaP,selection);
            phi_M = please_Give_Flux_Limiter(thetaM,selection);
            
            F(i,len) =  -lamY*cBz(len,i) - 0.5*lamY*(lamY-1)*( phi_P*cFz(len,i) - phi_M*cBz(len,i) );
        end

    end
    
    %Standard Godunov Flux Limiter
    for i=1:len
        for j=2:len-1
            if signs(i,j) <= 0 % u < 0
                
                if j==len-1
                    thetaP = ( C(1,i) - C(j+1,i) ) / ( C(j+1,i) - C(j,i) );
                    thetaM = ( C(j+1,i) - C(j,i) ) / ( C(j,i) - C(j-1,i) );                    
                else
                    thetaP = ( C(j+2,i) - C(j+1,i) ) / ( C(j+1,i) - C(j,i) );
                    thetaM = ( C(j+1,i) - C(j,i) ) / ( C(j,i) - C(j-1,i) );                   
                end
                
                lamY = uZ(j,i)*dt/dz;
                phi_P = please_Give_Flux_Limiter(thetaP,selection);
                phi_M = please_Give_Flux_Limiter(thetaM,selection);

                F(j,i) = -lamY*cFz(j,i) - 0.5*lamY*(lamY-1)*( phi_M*cBz(j,i) - phi_P*cFz(j,i) );
                
            else               % u > 0
                
                if j==2
                    thetaP = ( C(j,i) - C(j-1,i) ) / ( C(j+1,i) - C(j,i) );
                    thetaM = ( C(j-1,i) - C(len,i) ) / ( C(j,i) - C(j-1,i) );                      
                else
                    thetaP = ( C(j,i) - C(j-1,i) ) / ( C(j+1,i) - C(j,i) );
                    thetaM = ( C(j-1,i) - C(j-2,i) ) / ( C(j,i) - C(j-1,i) );   
                end
                
                lamY = uZ(j,i)*dt/dz;
                phi_P = please_Give_Flux_Limiter(thetaP,selection);
                phi_M = please_Give_Flux_Limiter(thetaM,selection);
                
                F(j,i) = -lamY*cFz(j,i) - 0.5*lamY*(lamY-1)*( phi_P*cFz(j,i) - phi_M*cBz(j,i) );
            end
        end
    end

    % Ends y-Direction calculation %
    
else
        
    fprintf('\n\n\n ERROR IN FUNCTION FOR COMPUTING FLUX LIMITERS\n');
       
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes flux limiter with choice of which one
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function phi = please_Give_Flux_Limiter(r,selection)

if strcmp(selection,'superbee')
    max1 = max( min(1,2*r),min(2,r) );
    phi = max(0,max1);
elseif strcmp(selection,'vanLeer')
    phi = ( r + abs(r) ) / ( 1 + abs(r) );
else
    fprintf('\n\n');
    error('NEED TO CHOOSE AN APPROPRIATE FLUX LIMITER');
end