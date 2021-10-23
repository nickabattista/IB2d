%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Setting up advection-diffusion solver
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C,laplacian_C] = please_Update_Adv_Diff_Concentration_Flux_Limiter_FV(C,dt,dx,dy,uX,uY,k)

% C:     concentration 
% dt:    time-step
% dx,dy: spatial steps in x and y, respectively
% uX:    x-Component of Velocity
% uY:    y-Component of Velocity
% k:     diffusion coefficient

% Compute Fluxes (Note: these calculations could be parallalized)
selection = 'superbee';
Fx = give_Necessary_Fluxes(C,dx,uX,'x',selection,dt); % Fluxes in x
Fy = give_Necessary_Fluxes(C,dy,uY,'y',selection,dt); % Fluxes in y

% "forward difference of fluxes in x"
F2x = [Fx(:,2:end) Fx(:,1)];
F1x = [Fx(:,end) Fx(:,1:end-1)];
diffX = 0.5/dx*( F2x - F1x );
   
% "forward differences of fluxes in y"
F2y = [Fy(2:end,:); Fy(1,:)];
F1y = [Fy(end,:); Fy(1:end-1,:)];
diffY = 0.5/dy*( F2y - F1y );
   

% Compute 2nd Derivative Terms
Cxx = DD(C,dx,'x');
Cyy = DD(C,dy,'y');


% Forms Laplacian
laplacian_C = Cxx+Cyy;
    
% UPWIND
C = C - dt * ( diffX + diffY ) + dt*( k*laplacian_C );



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Computes derivative based on sign of Velocity, u, using UPWIND
% approach
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C_z = give_Necessary_Fluxes(C,dz,uZ,string,selection,dt)

C_z = zeros(size(C));
len = length(uZ(:,1));
signs = sign(uZ);

if strcmp(string,'x')

    %For periodicity on ends w/ UPWIND
    for i=1:len

        %left side of grid
        if signs(i,1) <= 0 
            r = ( C(i,2) - C(i,1) ) / ( C(i,1) - C(i,len) );
            phi = please_Give_Flux_Limiter(r,selection);
            C_z(i,1) =  uZ(i,1)*C(i,2)   + 0.5*abs( uZ(i,1) )*( 1 - abs( 0.5*dt*uZ(i,1)/dz ) )*phi*( C(i,2) - C(i,len) );
        else
            r = ( C(i,1) - C(i,len) ) / ( C(i,2) - C(i,1) );
            phi = please_Give_Flux_Limiter(r,selection);
            C_z(i,1) =  uZ(i,1)*C(i,len) + 0.5*abs( uZ(i,1) )*( 1 - abs( 0.5*dt*uZ(i,1)/dz ) )*phi*( C(i,2) - C(i,len) );
        end

        %right side of grid
        if signs(len,1) <= 0
            r = ( C(i,1) - C(i,len) ) / ( C(i,len) - C(i,len-1) );
            phi = please_Give_Flux_Limiter(r,selection);
            C_z(i,len) =  uZ(i,len)*C(i,1)     + 0.5*abs( uZ(i,len) )*( 1 - 0.5*abs( dt*uZ(i,len)/dz ) )*phi*( C(i,1) - C(i,len-1) );
        else
            r = ( C(i,len) - C(i,len-1) ) / ( C(i,1) - C(i,len) );
            phi = please_Give_Flux_Limiter(r,selection);
            C_z(i,len) =  uZ(i,len)*C(i,len-1) + 0.5*abs( uZ(i,len) )*( 1 - 0.5*abs( dt*uZ(i,len)/dz ) )*phi*( C(i,1) - C(i,len-1) );
        end

    end
    %Standard Upwind 
    for i=1:len
        for j=2:len-1
            if signs(i,j) <= 0
                r = ( C(i,j+1) - C(i,j) ) / ( C(i,j) - C(i,j-1) );
                phi = please_Give_Flux_Limiter(r,selection);
                C_z(i,j) = uZ(i,j)*C(i,j+1) + 0.5*abs( uZ(i,j) )*( 1 - 0.5*abs( dt*uZ(i,j)/dz ) )*phi*( C(i,j+1) - C(i,j-1) );
            else
                r = ( C(i,j) - C(i,j-1) ) / ( C(i,j+1) - C(i,j) );
                phi = please_Give_Flux_Limiter(r,selection);
                C_z(i,j) = uZ(i,j)*C(i,j-1) + 0.5*abs( uZ(i,j) )*( 1 - 0.5*abs( dt*uZ(i,j)/dz ) )*phi*( C(i,j+1) - C(i,j-1) );
            end
        end
    end

    % Ends x-Direction calculation %
    
elseif strcmp(string,'y') 

    %For periodicity on ends w/ UPWIND
    for i=1:len

        %bottom of grid
        if signs(1,i) <= 0
            r = ( C(2,i) - C(1,i) ) / ( C(1,i) - C(len,i) );
            phi = please_Give_Flux_Limiter(r,selection);
            C_z(1,i) =  uZ(1,i)*C(2,i)   + 0.5*abs( uZ(1,i) )*( 1 - abs( 0.5*dt*uZ(1,i)/dz ) )*phi*( C(2,i) - C(len,i) );
        else
            r = ( C(1,i) - C(len,i) ) / ( C(2,i) - C(1,i) );
            phi = please_Give_Flux_Limiter(r,selection);
            C_z(1,i) =  uZ(1,i)*C(len,i) + 0.5*abs( uZ(1,i) )*( 1 - abs( 0.5*dt*uZ(1,i)/dz ) )*phi*( C(2,i) - C(len,i) );
        end

        %top of grid
        if signs(len,1) <= 0
            r = ( C(1,i) - C(len,i) ) / ( C(len,i) - C(len-1,i) );
            phi = please_Give_Flux_Limiter(r,selection);
            C_z(len,i) =  uZ(len,i)*C(1,i)       + 0.5*abs( uZ(len,i) )*( 1 - abs( 0.5*dt*uZ(len,i)/dz ) )*phi*( C(1,i) - C(len-1,i) );
        else
            r = ( C(len,i) - C(len-1,i) ) / ( C(1,i) - C(len,i) );
            phi = please_Give_Flux_Limiter(r,selection);
            C_z(len,i) =  uZ(len,i)*C(len-1,i)   + 0.5*abs( uZ(len,i) )*( 1 - abs( 0.5*dt*uZ(len,i)/dz ) )*phi*( C(1,i) - C(len-1,i) );
        end

    end
    
    %Standard Upwind
    for i=1:len
        for j=2:len-1
            if signs(j,i) <= 0
                r = ( C(j+1,i) - C(j,i) ) / ( C(j,i) - C(j-1,i) );
                phi = please_Give_Flux_Limiter(r,selection);
                C_z(j,i) = uZ(j,i)*C(j+1,i) + 0.5*abs( uZ(j,i) )*( 1 - abs( 0.5*dt*uZ(j,i)/dz ) )*phi*( C(j+1,i) - C(j-1,i) );
            else
                r = ( C(j,i) - C(j-1,i) ) / ( C(j+1,i) - C(j,i) );
                phi = please_Give_Flux_Limiter(r,selection);
                C_z(j,i) = uZ(j,i)*C(j-1,i) + 0.5*abs( uZ(j,i) )*( 1 - abs( 0.5*dt*uZ(j,i)/dz ) )*phi*( C(j+1,i) - C(j-1,i) );
            end
        end
    end

    % Ends y-Direction calculation %
    
else
        
    fprintf('\n\n\n ERROR IN FUNCTION FOR COMPUTING FLUX LIMITERS\n');
       
end
    
clear signs; clear len;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes flux limiter with choice of which one
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function phi = please_Give_Flux_Limiter(r,selection)

if strcmp(selection,'superbee')
    max1 = max( min(1,2*r),min(2,r) );
    phi = max(0,max1);
else
    fprintf('\n\n');
    error('NEED TO CHOOSE AN APPROPRIATE FLUX LIMITER');
end