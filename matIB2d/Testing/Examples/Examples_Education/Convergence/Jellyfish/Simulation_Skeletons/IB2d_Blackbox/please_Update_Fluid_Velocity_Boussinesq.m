%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  nick.battista@unc.edu
% Date Created: May 27th, 2015
% Institution: UNC-CH
%
% This code is capable of creating Lagrangian Structures using:
% 	1. Springs
% 	2. Beams (*torsional springs)
% 	3. Target Points
%	4. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-Tension)")
%
% One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting %%	lengths, etc
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
% 
% If you would like us %to add a specific muscle model, please let Nick (nick.battista@unc.edu) know.
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function solves the incompressible Navier-Stokes (NS) equations using Fast-Fourier Transform (FFT)
%      
%      x-Momentum Conservation: rho*u_t = -rho*u*u_x + rho*v*u_y + mu*laplacian(u) - p_x + Fx
%      y-Momentum Conservation: rho*v_t = -rho*u*v_x + rho*v*v_y + mu*laplacian(v) - p_y + Fy
%      Incompressibility:       u_x + v_y = 0.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [U_h, V_h, U, V, p] = please_Update_Fluid_Velocity_Boussinesq(U, V, Fx, Fy, rho, mu, grid_Info, dt)
 

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


% Construct EULERIAN Index Matrices
indy_X = (0:1:Nx-1);  idX = [];
indy_Y = (0:1:Ny-1)'; idY = [];
%Create x-Indices Grid
for i=1:Nx
    idX = [idX; indy_X]; 
end
%Create y-Indices Grid
for i=1:Ny
    idY = [idY indy_Y];
end


% Create FFT Operator (used for both half time-step and full time-step computations)
%A_hat = 1 + 2*mu*dt./rho'*( (sin(pi*idX/Nx)/dx).^2 + (sin(pi*idY/Ny)/dy).^2 );
A_hat = 1 + 2*mu*dt./rho*( (sin(pi*idX/Nx)/dx).^2 + (sin(pi*idY/Ny)/dy).^2 );


% % % % % % EVOLVE THE FLUID THROUGH HALF A TIME-STEP % % % % % % %


% Compute 1ST and 2ND derivatives (using centered differencing)
Ux = D(U,dx,'x');
Uy = D(U,dy,'y');

Uxx = DD(U,dx,'x');
Uyy = DD(U,dy,'y');

Vx = D(V,dx,'x');
Vy = D(V,dy,'y');

Vxx = DD(V,dx,'x');
Vyy = DD(V,dy,'y');

% Find derivatives of products U^2, V^2, and U*.
U_sq = U.^2;
V_sq = V.^2;
UV = U.*V;

U_sq_x = D(U_sq,dx,'x');
V_sq_y = D(V_sq,dy,'y');

UV_x = D(UV,dx,'x');
UV_y = D(UV,dy,'y');

%
% Construct right hand side in linear system
%
rhs_u = give_RHS_HALF_Step_Velocity(dt,rho,Nx,Ny,U,Ux,Uy,U_sq_x,UV_y,V,Fx,'x');
rhs_v = give_RHS_HALF_Step_Velocity(dt,rho,Nx,Ny,V,Vx,Vy,V_sq_y,UV_x,U,Fy,'y');


% Perform FFT to take velocities to state space
rhs_u_hat = fft2(rhs_u);
rhs_v_hat = fft2(rhs_v);  

% Calculate Fluid Pressure
p_hat = give_Fluid_Pressure(0.5*dt,rho,dx,dy,Nx,Ny,idX,idY,rhs_u_hat,rhs_v_hat);


% Calculate Fluid Velocity
u_hat = give_Me_Fluid_Velocity(0.5*dt,rho,dx,Nx,Ny,rhs_u_hat,p_hat,A_hat,idX,'x');
v_hat = give_Me_Fluid_Velocity(0.5*dt,rho,dy,Nx,Ny,rhs_v_hat,p_hat,A_hat,idY,'y');


% Inverse FFT to Get Velocities in Real Space
U_h = real(ifft2(u_hat));   %Half-step velocity, u
V_h = real(ifft2(v_hat));   %Half-step velocity, v




% % % % % % NOW EVOLVE THE FLUID THROUGH A FULL TIME-STEP % % % % % % %



% Compute first derivatives (centred) at half step.
U_h_x = D(U_h,dx,'x');
U_h_y = D(U_h,dy,'y');
V_h_x = D(V_h,dx,'x');
V_h_y = D(V_h,dy,'y');

% Computed derivatives of products U^2, V^2, and U*V at half step.
U_h_sq = U_h.^2;
V_h_sq = V_h.^2;
U_h_V_h = U_h.*V_h;
U_h_sq_x = D(U_h_sq,dx,'x');
V_h_sq_y = D(V_h_sq,dy,'y');
U_h_V_h_x = D(U_h_V_h,dx,'x');
U_h_V_h_y = D(U_h_V_h,dy,'y');

% Construct right hand side in linear system
%rhs_u = give_RHS_FULL_Step_Velocity(dt,mu,rho,Nx,Ny,U,Ux,Uy,U_sq_x,UV_y,V,Fx,Uxx,Uyy,'x');
%rhs_v = give_RHS_FULL_Step_Velocity(dt,mu,rho,Nx,Ny,V,Vx,Vy,V_sq_y,UV_x,U,Fy,Vxx,Vyy,'y');
rhs_u = give_RHS_FULL_Step_Velocity(dt,mu,rho,Nx,Ny,U,U_h_x,U_h_y,U_h_sq_x,U_h_V_h_y,V,Fx,Uxx,Uyy,'x');
rhs_v = give_RHS_FULL_Step_Velocity(dt,mu,rho,Nx,Ny,V,V_h_x,V_h_y,V_h_sq_y,U_h_V_h_x,U,Fy,Vxx,Vyy,'y');


% Perform FFT to take velocities to state space
rhs_u_hat = fft2(rhs_u);
rhs_v_hat = fft2(rhs_v);  


% Calculate Fluid Pressure
p_hat  = give_Fluid_Pressure(dt,rho,dx,dy,Nx,Ny,idX,idY,rhs_u_hat,rhs_v_hat);
        

% Calculate Fluid Velocity
u_hat = give_Me_Fluid_Velocity(dt,rho,dx,Nx,Ny,rhs_u_hat,p_hat,A_hat,idX,'x');
v_hat = give_Me_Fluid_Velocity(dt,rho,dy,Nx,Ny,rhs_v_hat,p_hat,A_hat,idY,'y');


% Inverse FFT to Get Velocities/Pressure in Real Space
U = real(ifft2(u_hat));
V = real(ifft2(v_hat));
p = real(ifft2(p_hat));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: calculates the fluid velocity!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vel_hat = give_Me_Fluid_Velocity(dt,rho,dj,Nx,Ny,rhs_VEL_hat,p_hat,A_hat,idMat,string)

vel_hat = zeros(Ny,Nx); %initialize fluid velocity


if strcmp(string,'x')
    
    for i=1:Ny
        for j=1:Nx
            vel_hat(i,j) = ( rhs_VEL_hat(i,j) - 1i*dt/(rho(i,j)*dj)*sin(2*pi*idMat(i,j)/Nx)*p_hat(i,j) ) / A_hat(i,j);
            %vel_hat(i,j) = ( rhs_VEL_hat(i,j) - 1i*dt/(rho(j,i)*dj)*sin(2*pi*idMat(i,j)/Nx)*p_hat(i,j) ) / A_hat(i,j);

        end
    end

elseif strcmp(string,'y')
    
    for i=1:Ny
        for j=1:Nx
            vel_hat(i,j) = ( rhs_VEL_hat(i,j) - 1i*dt./(rho(i,j)*dj)*sin(2*pi*idMat(i,j)/Ny)*p_hat(i,j) ) / A_hat(i,j);
            %vel_hat(i,j) = ( rhs_VEL_hat(i,j) - 1i*dt./(rho(j,i)*dj)*sin(2*pi*idMat(i,j)/Ny)*p_hat(i,j) ) / A_hat(i,j);

        end
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates RHS with fluid velocity in FULL STEP computation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rhs = give_RHS_FULL_Step_Velocity(dt,mu,rho,Nx,Ny,A,Ax,Ay,A_sq_j,AB_j,B,Fj,Axx,Ayy,string)

% Note: Fj -> j corresponds to either x or y.

rhs = zeros(Ny,Nx); %initialize rhs

if strcmp(string,'x')
    
    for i=1:Ny
        for j=1:Nx
            rhs(i,j) = A(i,j) + dt/rho(i,j)*( Fj(i,j) + mu/2*(Axx(i,j)+Ayy(i,j)) - 0.5*rho(i,j)*(A(i,j)*Ax(i,j) + B(i,j)*Ay(i,j)) - .5*rho(i,j)*(A_sq_j(i,j) + AB_j(i,j) ) ); %RHS: u-component
            %rhs(i,j) = A(i,j) + dt/rho(j,i)*( Fj(i,j) + mu/2*(Axx(i,j)+Ayy(i,j)) - 0.5*rho(j,i)*(A(i,j)*Ax(i,j) + B(i,j)*Ay(i,j)) - .5*rho(j,i)*(A_sq_j(i,j) + AB_j(i,j) ) ); %RHS: u-component

        end
    end
    
elseif strcmp(string,'y')
    
    for i=1:Ny
        for j=1:Nx
            rhs(i,j) = A(i,j) + dt/rho(i,j)*( Fj(i,j) + mu/2*(Axx(i,j)+Ayy(i,j)) - 0.5*rho(i,j)*(B(i,j)*Ax(i,j) + A(i,j)*Ay(i,j)) - .5*rho(i,j)*(AB_j(i,j) + A_sq_j(i,j) ) ); %RHS: v-compoent
            %rhs(i,j) = A(i,j) + dt/rho(j,i)*( Fj(i,j) + mu/2*(Axx(i,j)+Ayy(i,j)) - 0.5*rho(j,i)*(B(i,j)*Ax(i,j) + A(i,j)*Ay(i,j)) - .5*rho(j,i)*(AB_j(i,j) + A_sq_j(i,j) ) ); %RHS: v-compoent
        end
    end
   
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates RHS with fluid velocity in HALF computation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rhs = give_RHS_HALF_Step_Velocity(dt,rho,Nx,Ny,A,Ax,Ay,A_sq_j,AB_j,B,Fj,string)

% Note: Fj -> j corresponds to either x or y.

rhs = zeros(Ny,Nx); %initialize rhs


if strcmp(string,'x')
    
    for i=1:Ny
        for j=1:Nx
            rhs(i,j) = A(i,j) + 0.5*dt/rho(i,j)*( Fj(i,j) - 0.5*rho(i,j)*(A(i,j)*Ax(i,j) + B(i,j)*Ay(i,j)) - 0.5*rho(i,j)*(A_sq_j(i,j) + AB_j(i,j) ) ); %RHS: u-component
            %rhs(i,j) = A(i,j) + 0.5*dt/rho(j,i)*( Fj(i,j) - 0.5*rho(j,i)*(A(i,j)*Ax(i,j) + B(i,j)*Ay(i,j)) - 0.5*rho(j,i)*(A_sq_j(i,j) + AB_j(i,j) ) ); %RHS: u-component       
        end
    end
    
elseif strcmp(string,'y')
    
    for i=1:Ny
        for j=1:Nx
            rhs(i,j) = A(i,j) + 0.5*dt/rho(i,j)*( Fj(i,j) - 0.5*rho(i,j)*(B(i,j)*Ax(i,j) + A(i,j)*Ay(i,j)) - 0.5*rho(i,j)*(AB_j(i,j) + A_sq_j(i,j) ) ); %RHS: v-compoent
            %rhs(i,j) = A(i,j) + 0.5*dt/rho(j,i)*( Fj(i,j) - 0.5*rho(j,i)*(B(i,j)*Ax(i,j) + A(i,j)*Ay(i,j)) - 0.5*rho(j,i)*(AB_j(i,j) + A_sq_j(i,j) ) ); %RHS: v-compoent
        end
    end
   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: calculates the fluid pressure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p_hat = give_Fluid_Pressure(dt,rho,dx,dy,Nx,Ny,idX,idY,rhs_u_hat,rhs_v_hat)

p_hat = zeros(Ny,Nx); %initialize fluid pressure

for i=1:Ny
    for j=1:Nx
        num = -( 1i/dx*sin(2*pi*idX(i,j)/Nx)*rhs_u_hat(i,j) + 1i/dy*sin(2*pi*idY(i,j)/Ny)*rhs_v_hat(i,j) );
        den = ( dt/rho(i,j)*( (sin(2*pi*idX(i,j)/Nx)/dx)^2 + (sin(2*pi*idY(i,j)/Ny)/dy).^2 ));
        %den = ( dt/rho(j,i)*( (sin(2*pi*idX(i,j)/Nx)/dx)^2 + (sin(2*pi*idY(i,j)/Ny)/dy).^2 ));
        p_hat(i,j) = num/den;
    end 
end

% Zero out modes.
p_hat(1,1) = 0;
p_hat(1,Nx/2+1) = 0;
p_hat(Ny/2+1,Nx/2+1) = 0;  
p_hat(Ny/2+1,1) = 0;



