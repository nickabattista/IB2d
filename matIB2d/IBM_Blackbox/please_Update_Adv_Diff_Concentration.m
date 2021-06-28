%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Sets up the advection/diffusion solver for a background
%           concentration 
%
%   Upwind w/ Time Splitting
%   Author: Nick Battista
%   Date: June 2016 (created); June 2021 (modified)
%   Institution (when created):  UNC-CH
%   Current Instituteion: TCNJ
%
%   WENO SCHEMES
%   Author: Matea Santiago
%   Date: February 2021
%   Institution (when created): UC Merced
%
%   Returns: spatial derivatives in x and y of concentration
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C,laplacian_C] = please_Update_Adv_Diff_Concentration(C,dt,dx,dy,uX,uY,k,Adv_flag,Lx,Ly)

% C:     concentration 
% dt:    time-step
% dx,dy: spatial steps in x and y, respectively
% uX:    x-Component of Velocity
% uY:    y-Component of Velocity
% k:     diffusion coefficient


%
% Choose whether UPWIND or WENO advection is used:
%
if Adv_flag == 0
    
    % Performs Upwind Advection WITH Time-Splitting
    C = perform_Time_Split_Upwind(C,dt,dx,dy,uX,uY,k);
        
    % Performs Upwind Advection WITHOUT Time-Splitting
    %C = perform_Time_noSplit_Upwind(C,dt,dx,dy,uX,uY,k);

elseif Adv_flag == 1

    % Performs WENO Advection WITHOUT Time-Splitting
    C = perform_WENO_source(C,dt,dx,dy,uX,uY,k,0,Lx,Ly);

else
    'UNRECOGNIZED ADVECTION SCHEME'
end

% Performs Upwind Advection w/ Time-Splitting

laplacian_C=1; % DUMMY VARIABLE (laplacian not used anywhere else in code.)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Advection-Diffusion Split Upwind Method
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C = perform_Time_noSplit_Upwind(C,dt,dx,dy,uX,uY,k)

% Compute Necessary Derivatives (Note: these calculations could be parallalized)
Cx = give_Necessary_Derivative(C,dx,uX,'x');
Cy = give_Necessary_Derivative(C,dy,uY,'y'); 
Cxx = DD(C,dx,'x');
Cyy = DD(C,dy,'y');

% Forms Laplacian
laplacian_C = Cxx+Cyy;
    
% UPWIND
C = C + dt * ( k*(laplacian_C) - uX.*Cx - uY.*Cy );



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Advection-Diffusion Split Upwind Method
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C = perform_Time_Split_Upwind(C,dt,dx,dy,uX,uY,k)

% Compute Necessary Derivatives for x-Advection 
Cx = give_Necessary_Derivative(C,dx,uX,'x');
Cxx = DD(C,dx,'x');

% Time-step #1 (give auxillary)
C = C + dt * ( k*(Cxx) - uX.*Cx  );

% Compute Necessary Derivatives for y-Advection 
Cy = give_Necessary_Derivative(C,dy,uY,'y'); 
Cyy = DD(C,dy,'y');
 
% Time-step #2 (give next iteration)
C = C + dt * ( k*(Cyy) - uY.*Cy );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Advection-Diffusion w/ 3RD-order WENO Scheme for advection
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C = perform_WENO_source(C,dt,dx,dy,uX,uY,k,f,Lx,Ly)

% Compute Necessary Derivatives using WENO
[Cx,Cy]=WENO_3O(C,uX,uY,dx,dy,dt,Lx,Ly);


Cxx = DD(C,dx,'x');
Cyy = DD(C,dy,'y');

% Forms Laplacian
laplacian_C = Cxx+Cyy;
    
% Update
C=C+dt*(k*laplacian_C-(uX.*Cx+uY.*Cy)+f);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Computes derivative based on sign of Velocity, u, using UPWIND
% approach
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C_z = give_Necessary_Derivative(C,dz,uZ,string)

C_z = zeros(size(C));
signs = sign(uZ);
[lenY,lenX] = size(uZ);

if strcmp(string,'x')

    %For periodicity on ends w/ UPWIND
    for i=1:lenY

        %left side of grid
        if signs(i,1) <= 0 
            C_z(i,1) =  ( C(i,2) - C(i,1) ) / (dz);
        else
            C_z(i,1) =  ( C(i,1) - C(i,lenX) ) / (dz);
        end

        %right side of grid
        if signs(i,lenX) <= 0 
            C_z(i,lenX) =  ( C(i,1) - C(i,lenX) ) / (dz);
        else
            C_z(i,lenX) =  ( C(i,lenX) - C(i,lenX-1) ) / (dz);
        end

    end
    %Standard Upwind 
    for i=1:lenY
        for j=2:lenX-1
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
    for i=1:lenX

        %bottom of grid
        if signs(1,i) <= 0 
            C_z(1,i) =  ( C(2,i) - C(1,i) ) / (dz);
        else
            C_z(1,i) =  ( C(1,i) - C(lenY,i) ) / (dz);
        end

        %top of grid
        if signs(lenY,i) <= 0 
            C_z(lenY,i) =  ( C(1,i) - C(lenY,i) ) / (dz);
        else
            C_z(lenY,i) =  ( C(lenY,i) - C(lenY-1,i) ) / (dz);
        end

    end
    
    %Standard Upwind
    for i=2:lenY-1
        for j=1:lenX
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

