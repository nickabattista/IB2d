%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  battistn[@]tcnj[.]edu
% 
% IB2d was Created: May 27th, 2015 at UNC-CH
%
% This code is capable of creating Lagrangian Structures using:
% 	1. Springs
% 	2. Beams (*torsional springs*)
% 	3. Target Points
% 	4. Muscle-Model (combined Force-Length-Velocity model, "Hill+(Length-Tension)")
%   .
%   .
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function solves the incompressible Navier-Stokes (NS) equations using Fast-Fourier Transform (FFT)
%      
%      x-Momentum Conservation: rho*u_t = -rho*u*u_x + rho*v*u_y + mu*laplacian(u) - p_x + Fx
%      y-Momentum Conservation: rho*v_t = -rho*u*v_x + rho*v*v_y + mu*laplacian(v) - p_y + Fy
%      Incompressibility:       u_x + v_y = 0
%
%      NOTE: (i) Lots of old implementation included (but commented out) for teaching purposes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function [U_h, V_h, U, V, p] = please_Update_Fluid_Velocity(U, V, Fx, Fy, rho, mu, grid_Info, dt, idX, idY, A_hat)
function [U_h, V_h, U, V, p] = please_Update_Fluid_Velocity(U, V, Fx, Fy, rho, mu, grid_Info, dt, sinIDX, sinIDY, A_hat)
 

% Fluid (Eulerian) Grid updated using Peskin's two-step algorithm, where the advection terms
% are written in skew symmetric form.
%
% U:         Eulerian grid of x-Velocities
% V:         Eulerian grid of y-Velocities
% Fx:        x-Forces arising from deformations in Lagrangian structure
% Fy:        y-Forces arising from deformatinos in Lagrangian structure
% rho:       Fluid density
% mu:        Fluid dynamic viscosity
% grid_Info: Vector of parameters relating to Eulerian grid, Lagrangian grid, numerical delta function
% dt:        Time-step
% idX/idY:   EULERIAN Index Matrices for FFT Operators


% Initialize %
Nx =   grid_Info(1);
Ny =   grid_Info(2);
Lx =   grid_Info(3);
Ly =   grid_Info(4);
dx =   grid_Info(5);
dy =   grid_Info(6);
supp = grid_Info(7);
Nb =   grid_Info(8);
ds =   grid_Info(9);



%-----------------------------------------------------------------------------
%           **** NOW INITIALIZED BEFORE TIME-LOOP STARTS!!! ****
% Create FFT Operator (used for both half time-step and full time-step computations)
%-----------------------------------------------------------------------------
%A_hat = 1 + 2*mu*dt/rho*( (sin(pi*idX/Nx)/dx).^2 + (sin(pi*idY/Ny)/dy).^2 );


%**************************************************************************
%**************************************************************************
%
% % % % % % EVOLVE THE FLUID THROUGH HALF A TIME-STEP % % % % % % %
%
%**************************************************************************
%**************************************************************************


%----------------------------------------------------------------
%        MORE EFFICIENT IMPLEMENTATIONS OF D() and DD()
% Compute 1ST and 2ND derivatives (using centered differencing)
%----------------------------------------------------------------
Ux = D(U,dx,'x');
Uy = D(U,dy,'y');

Uxx = DD(U,dx,'x');
Uyy = DD(U,dy,'y');

Vx = D(V,dx,'x');
Vy = D(V,dy,'y');

Vxx = DD(V,dx,'x');
Vyy = DD(V,dy,'y');

%-----------------------------------------------------
%         ORIGINAL SLOWER IMPLEMENTATION!!!
% Find derivatives of products U^2, V^2, and U.*V
%-----------------------------------------------------
% U_sq = U.^2;
% V_sq = V.^2;
% UV = U.*V;
% U_sq_x2 = D(U_sq,dx,'x');
% V_sq_y2 = D(V_sq,dy,'y');
% UV_x2 = D(UV,dx,'x');
% UV_y2 = D(UV,dy,'y');

%----------------------------------------------------------------
%               MORE EFFICIENT IMPLEMENTATION
% Find derivatives of products U^2, V^2, and U.*V at half step
%               using product and chain rules
%----------------------------------------------------------------
U_sq_x = 2*U.*Ux;
V_sq_y = 2*V.*Vx;
%
UV_x = V.*Ux + U.*Vx;
UV_y = V.*Uy + U.*Vy;


%-----------------------------------------------------
% Construct right hand side in linear system
%-----------------------------------------------------
rhs_u = give_RHS_HALF_Step_Velocity(dt,rho,Nx,Ny,U,Ux,Uy,U_sq_x,UV_y,V,Fx,'x');
rhs_v = give_RHS_HALF_Step_Velocity(dt,rho,Nx,Ny,V,Vx,Vy,V_sq_y,UV_x,U,Fy,'y');


%-----------------------------------------------------
% Perform FFT to take velocities to state space
%-----------------------------------------------------
rhs_u_hat = fft2(rhs_u);
rhs_v_hat = fft2(rhs_v);  


%--------------------------------------------------------------
% Calculate Fluid Pressure (uses stored FOURIER matrices)
%--------------------------------------------------------------
p_hat = give_Fluid_Pressure(0.5*dt,rho,dx,dy,Nx,Ny,sinIDX,sinIDY,rhs_u_hat,rhs_v_hat);


%--------------------------------------------------------------
% Calculate Fluid Velocity (USING stored SINE vals)
%--------------------------------------------------------------
u_hat = give_Me_Fluid_Velocity(0.5*dt,rho,dx,Nx,Ny,rhs_u_hat,p_hat,A_hat,sinIDX,'x');
v_hat = give_Me_Fluid_Velocity(0.5*dt,rho,dy,Nx,Ny,rhs_v_hat,p_hat,A_hat,sinIDY,'y');


%--------------------------------------------------------------
% Inverse FFT to Get Velocities in Real Space
%--------------------------------------------------------------
U_h = real(ifft2(u_hat));   % Half-step velocity, u
V_h = real(ifft2(v_hat));   % Half-step velocity, v


%**************************************************************************
%**************************************************************************
%
% % % % % % NOW EVOLVE THE FLUID THROUGH A FULL TIME-STEP % % % % % % %
%
%**************************************************************************
%**************************************************************************

%-----------------------------------------------------------------
%        MORE EFFICIENT IMPLEMENTATIONS OF D() and DD()
%      Compute first derivatives (centered) at half step.
%-----------------------------------------------------------------
U_h_x = D(U_h,dx,'x');
U_h_y = D(U_h,dy,'y');
V_h_x = D(V_h,dx,'x');
V_h_y = D(V_h,dy,'y');


%-----------------------------------------------------------------
%              ORIGINAL SLOWER IMPLEMENTATION!!!
% Find derivatives of products U^2, V^2, and U.*V at half step.
%-----------------------------------------------------------------
% U_h_sq = U_h.^2;
% V_h_sq = V_h.^2;
% U_h_V_h = U_h.*V_h;
% U_h_sq_x2 = D(U_h_sq,dx,'x');
% V_h_sq_y2 = D(V_h_sq,dy,'y');
% U_h_V_h_x2 = D(U_h_V_h,dx,'x');
% U_h_V_h_y2 = D(U_h_V_h,dy,'y');


%----------------------------------------------------------------
%               MORE EFFICIENT IMPLEMENTATION
% Find derivatives of products U^2, V^2, and U.*V at half step
%               using product and chain rules
%----------------------------------------------------------------
U_h_sq_x = 2*U_h.*U_h_x;
V_h_sq_y = 2*V_h.*V_h_y;
U_h_V_h_x = (U_h_x .* V_h) + (U_h .* V_h_x);
U_h_V_h_y = (U_h_y .* V_h) + (U_h .* V_h_y);


%---------------------------------------------------
% Construct right hand side in linear system
%---------------------------------------------------
%rhs_u = give_RHS_FULL_Step_Velocity(dt,mu,rho,Nx,Ny,U,Ux,Uy,U_sq_x,UV_y,V,Fx,Uxx,Uyy,'x');
%rhs_v = give_RHS_FULL_Step_Velocity(dt,mu,rho,Nx,Ny,V,Vx,Vy,V_sq_y,UV_x,U,Fy,Vxx,Vyy,'y');
rhs_u = give_RHS_FULL_Step_Velocity(dt,mu,rho,Nx,Ny,U,U_h_x,U_h_y,U_h_sq_x,U_h_V_h_y,V,Fx,Uxx,Uyy,'x');
rhs_v = give_RHS_FULL_Step_Velocity(dt,mu,rho,Nx,Ny,V,V_h_x,V_h_y,V_h_sq_y,U_h_V_h_x,U,Fy,Vxx,Vyy,'y');


%---------------------------------------------------
% Perform FFT to take velocities to state space
%---------------------------------------------------
rhs_u_hat = fft2(rhs_u);
rhs_v_hat = fft2(rhs_v);  


%-------------------------------------------------------------
% Calculate Fluid Pressure (uses stored FOURIER matrices)
%-------------------------------------------------------------
p_hat  = give_Fluid_Pressure(dt,rho,dx,dy,Nx,Ny,sinIDX,sinIDY,rhs_u_hat,rhs_v_hat);
        

%-------------------------------------------------------------
% Calculate Fluid Velocity (using stored sine values)
%-------------------------------------------------------------
u_hat = give_Me_Fluid_Velocity(dt,rho,dx,Nx,Ny,rhs_u_hat,p_hat,A_hat,sinIDX,'x');
v_hat = give_Me_Fluid_Velocity(dt,rho,dy,Nx,Ny,rhs_v_hat,p_hat,A_hat,sinIDY,'y');

%-------------------------------------------------------------
% Inverse FFT to Get Velocities/Pressure in Real Space
%-------------------------------------------------------------
U = real(ifft2(u_hat));
V = real(ifft2(v_hat));
p = real(ifft2(p_hat));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: calculates the fluid velocity!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function vel_hat = give_Me_Fluid_Velocity(dt,rho,dj,Nx,Ny,rhs_VEL_hat,p_hat,A_hat,idMat,string)
function vel_hat = give_Me_Fluid_Velocity(dt,rho,dj,Nx,Ny,rhs_VEL_hat,p_hat,A_hat,sinIDXorIDY,string)


%----------------------------------------------------------------------
% No initialization necessary since using vectorized operations 
%----------------------------------------------------------------------
%vel_hat = zeros(Ny,Nx); %initialize fluid velocity

if strcmp(string,'x')
    
      % SLOWER (for-loop, non-vectorized computation)
      %for i=1:Ny
      %    for j=1:Nx
      %        vel_hat(i,j) = ( rhs_VEL_hat(i,j) - 1i*dt/(rho*dj)*sin(2*pi*idMat(i,j)/Nx)*p_hat(i,j) ) / A_hat(i,j);
      %    end
      %end
    
      % VECTORIZED function for speedup (not using stored FOURIER/SINE VALS):
      %vel_hat = ( rhs_VEL_hat - 1i*dt/(rho*dj)*sin(2*pi*idMat/Nx).*p_hat ) ./ A_hat;
      
      % VECTORIZED w/ STORED SINE VALS
      vel_hat = ( rhs_VEL_hat - 1i*dt/(rho*dj)* sinIDXorIDY .*p_hat ) ./ A_hat;

elseif strcmp(string,'y')
    
      % SLOWER (for-loop, non-vectorized computation)
      %for i=1:Ny
      %    for j=1:Nx
      %        vel_hat(i,j) = ( rhs_VEL_hat(i,j) - 1i*dt/(rho*dj)*sin(2*pi*idMat(i,j)/Ny)*p_hat(i,j) ) / A_hat(i,j);
      %    end
      %end

      % VECTORIZED function for speedup (not using stored FOURIER/SINE VALS):
      %vel_hat = ( rhs_VEL_hat - 1i*dt/(rho*dj)*sin(2*pi*idMat/Ny).*p_hat ) ./ A_hat;
      
      % VECTORIZED w/ STORED SINE VALS
      vel_hat = ( rhs_VEL_hat - 1i*dt/(rho*dj)* sinIDXorIDY .*p_hat ) ./ A_hat;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates RHS with fluid velocity in FULL STEP computation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rhs = give_RHS_FULL_Step_Velocity(dt,mu,rho,Nx,Ny,A,Ax,Ay,A_sq_j,AB_j,B,Fj,Axx,Ayy,string)

% Note: Fj -> j corresponds to either x or y.

%----------------------------------------------------------------------
% No initialization necessary since using vectorized operations 
%----------------------------------------------------------------------
%rhs = zeros(Ny,Nx); %initialize rhs

if strcmp(string,'x')
    
    % SLOWER (for-loop, non-vectorized computation)
    %for i=1:Ny
    %    for j=1:Nx
    %       rhs(i,j) = A(i,j) + dt/rho*( Fj(i,j) + mu/2*(Axx(i,j)+Ayy(i,j)) - 0.5*rho*(A(i,j)*Ax(i,j) + B(i,j)*Ay(i,j)) - .5*rho*(A_sq_j(i,j) + AB_j(i,j) ) ); %RHS: u-component
    %    end
    %end
    
    % VECTORIZED function for speedup:
    rhs = A + (dt/rho)*( Fj + mu/2*(Axx+Ayy) - 0.5*rho*(A.*Ax + B.*Ay) - 0.5*rho*(A_sq_j + AB_j ) ); %RHS: u-component
    
elseif strcmp(string,'y')
    
    % SLOWER (for-loop, non-vectorized computation)
    %for i=1:Ny
    %    for j=1:Nx
    %        rhs(i,j) = A(i,j) + dt/rho*( Fj(i,j) + mu/2*(Axx(i,j)+Ayy(i,j)) - 0.5*rho*(B(i,j)*Ax(i,j) + A(i,j)*Ay(i,j)) - .5*rho*(AB_j(i,j) + A_sq_j(i,j) ) ); %RHS: v-component
    %    end
    %end
    
    % VECTORIZED function for speedup:
    rhs = A + (dt/rho)*( Fj + mu/2*(Axx+Ayy) - 0.5*rho*(B.*Ax + A.*Ay) - 0.5*rho*(AB_j + A_sq_j ) ); %RHS: v-component
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates RHS with fluid velocity in HALF computation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rhs = give_RHS_HALF_Step_Velocity(dt,rho,Nx,Ny,A,Ax,Ay,A_sq_j,AB_j,B,Fj,string)

% Note: Fj -> j corresponds to either x or y.

%----------------------------------------------------------------------
% No initialization necessary since using vectorized operations 
%----------------------------------------------------------------------
%rhs = zeros(Ny,Nx); %initialize rhs

if strcmp(string,'x')
    
    % SLOWER (for-loop, non-vectorized computation)
    %for i=1:Ny
    %    for j=1:Nx
    %        rhs(i,j) = A(i,j) + .5*dt/rho*( Fj(i,j) - .5*rho*(A(i,j)*Ax(i,j) + B(i,j)*Ay(i,j)) - .5*rho*(A_sq_j(i,j) + AB_j(i,j) ) ); %RHS: u-component
    %    end
    %end
    
    % VECTORIZED function for speedup:
    rhs = A + (0.5*dt/rho)*( Fj - 0.5*rho*(A.*Ax + B.*Ay) - 0.5*rho*( A_sq_j + AB_j ) ); %RHS: u-component
    
elseif strcmp(string,'y')
    
    % SLOWER (for-loop, non-vectorized computation)
    %for i=1:Ny
    %    for j=1:Nx
    %        rhs(i,j) = A(i,j) + .5*dt/rho*( Fj(i,j) - .5*rho*(B(i,j)*Ax(i,j) + A(i,j)*Ay(i,j)) - .5*rho*(AB_j(i,j) + A_sq_j(i,j) ) ); %RHS: v-compoent
    %    end
    %end
    
    % VECTORIZED function for speedup:
    rhs = A + (0.5*dt/rho)*( Fj - 0.5*rho*( B.*Ax + A.*Ay ) - .5*rho*( AB_j + A_sq_j ) ); %RHS: v-compoent
   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: calculates the fluid pressure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function p_hat = give_Fluid_Pressure(dt,rho,dx,dy,Nx,Ny,idX,idY,rhs_u_hat,rhs_v_hat)
function p_hat = give_Fluid_Pressure(dt,rho,dx,dy,Nx,Ny,sinIDX,sinIDY,rhs_u_hat,rhs_v_hat)

%----------------------------------------------------------------------
% No initialization necessary since using vectorized operations 
%----------------------------------------------------------------------
%p_hat = zeros(Ny,Nx); %initialize fluid pressure


% SLOWER (for-loop, non-vectorized computation)
%for i=1:Ny
%    for j=1:Nx
%        num = -( 1i/dx*sin(2*pi*idX(i,j)/Nx)*rhs_u_hat(i,j) + 1i/dy*sin(2*pi*idY(i,j)/Ny)*rhs_v_hat(i,j) );
%        den = ( dt/rho*( (sin(2*pi*idX(i,j)/Nx)/dx)^2 + (sin(2*pi*idY(i,j)/Ny)/dy).^2 ));
%        p_hat(i,j) = num/den;
%    end 
%end

% VECTORIZED function calculations for speedup:
%num = -( 1i/dx*sin(2*pi*idX/Nx).*rhs_u_hat + 1i/dy*sin(2*pi*idY/Ny).*rhs_v_hat );
%den = ( dt/rho*( ( sin(2*pi*idX/Nx)/dx ).^2 + ( sin(2*pi*idY/Ny)/dy ).^2 ) );

% ORIGINAL (not using stored SINE/FOURIER values)
%p_hat = ( -( 1i/dx*sin(2*pi*idX/Nx).*rhs_u_hat + 1i/dy*sin(2*pi*idY/Ny).*rhs_v_hat ) )./ ( dt/rho*( ( sin(2*pi*idX/Nx)/dx ).^2 + ( sin(2*pi*idY/Ny)/dy ).^2 ) );

% WITH STORED SINES
p_hat = ( -( 1i/dx* sinIDX .*rhs_u_hat + 1i/dy* sinIDY .*rhs_v_hat ) )./ ( dt/rho*( ( sinIDX /dx ).^2 + ( sinIDY /dy ).^2 ) );

% Zero out modes.
p_hat(1,1) = 0;
p_hat(1,Nx/2+1) = 0;
p_hat(Ny/2+1,Nx/2+1) = 0;  
p_hat(Ny/2+1,1) = 0;
