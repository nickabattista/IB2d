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
%   4. Mass Points
%   5. Porous Points
%	6. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-Tension)")
%   7. 3-Element Hill Muscle Model
%
% One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting-lengths, etc
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
% 
% If you would like us to add a specific muscle model, please let Nick (nick.battista@unc.edu) know.
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Actual DRIVER of the code, where the time-stepping occurs ->
%           gets called by main2d to do the "magic" :)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X, Y, U, V, xLag, yLag] = IBM_Driver(struct_name, mu, rho, grid_Info, dt, T_FINAL, model_Info)

%
%    2D IMMERSED BOUNDARY SOLVER ON RECTANGULAR DOMAIN w/ PERIODIC BOUNDARIES
%    
%    x-Momentum Conservation: rho*u_t = -rho*u*u_x + rho*v*u_y + mu*laplacian(u) - p_x + F_x
%    y-Momentum Convervation: rho*v_t = -rho*u*v_x + rho*v*v_y + mu*laplacian(v) - p_y + F_y
%
%    Incompressibility: u_x + v_y = 0
%
%    LagPts/dt = int{ u(x,t) delta( x - LagPts(s,t) ) dx }
%    F_x = int{ fx(s,t) delta(x - LagPts(s,t)) ds }
%    F_y = int{ fy(s,t) delta(x - LagPts(s,t)) ds }
%
fprintf('\n\n |****** Prepping Immersed Boundary Simulation ******|\n');
fprintf('\n\n--> Reading input data for simulation...\n\n');


% Temporal Information
NTime = floor(T_FINAL/dt)+1; % # of total time-steps (floor'd so exact number of time-steps)
dt = T_FINAL/NTime;          %time-step (slightly perturbed dt, so exact # of time-steps are used)
current_time = 0.0;          %initialize start of simulation to time, 0 


% GRID INFO %
Nx =   grid_Info(1); % # of Eulerian pts. in x-direction
Ny =   grid_Info(2); % # of Eulerian pts. in y-direction
Lx =   grid_Info(3); % Length of Eulerian grid in x-coordinate
Ly =   grid_Info(4); % Length of Eulerian grid in y-coordinate
dx =   grid_Info(5); % Spatial-size in x
dy =   grid_Info(6); % Spatial-size in y
supp = grid_Info(7); % Delta-function support

% PRINTING/PLOTTING INFO %
pDump = grid_Info(8);           % Print (Plot) Dump interval
pMatlab = grid_Info(9);         % Plot in Matlab? (1=YES,0=NO)
lagPlot = grid_Info(10);        % Plot LAGRANGIAN PTs ONLY in Matlab
velPlot = grid_Info(11);        % Plot LAGRANGIAN PTs + VELOCITY FIELD in Matlab
vortPlot = grid_Info(12);       % Plot LAGRANGIAN PTs + VORTICITY colormap in Matlab
uMagPlot = grid_Info(13);       % Plot LAGRANGIAN PTs + MAGNITUDE OF VELOCITY colormap in Matlab
pressPlot = grid_Info(14);      % Plot LAGRANGIAN PTs + PRESSURE colormap in Matlab


% MODEL STRUCTURE DATA STORED %
springs_Yes = model_Info(1);           % Springs: 0 (for no) or 1 (for yes) 
update_Springs_Flag = model_Info(2);   % Update_Springs: 0 (for no) or 1 (for yes)
target_pts_Yes = model_Info(3);        % Target_Pts: 0 (for no) or 1 (for yes)
update_Target_Pts = model_Info(4);     % Update_Target_Pts: 0 (for no) or 1 (for yes)
beams_Yes = model_Info(5);             % Beams: 0 (for no) or 1 (for yes)
update_Beams_Flag = model_Info(6);     % Update_Beams: 0 (for no) or 1 (for yes)
muscles_Yes = model_Info(7);           % FV-LT Muscles: 0 (for no) or 1 (for yes)
hill_3_muscles_Yes = model_Info(8);    % Hill 3-Element Muscle: 0 (for no) or 1 (for yes)
arb_ext_force_Yes = model_Info(9);     % Arbitrary External Force: 0 (for no) or 1 (for yes)
tracers_Yes = model_Info(10);          % Tracers: 0 (for no) or 1 (for yes)
mass_Yes = model_Info(11);             % Mass Points: 0 (for no) or 1 (for yes)
gravity_Yes = model_Info(12);          % Gravity: 0 (for no) or 1 (for yes)
%NOTE: model_Info(13)/(14) - components of gravity vector
porous_Yes = model_Info(15);           % Porous Media: 0 (for no) or 1 (for yes)
concentration_Yes = model_Info(16);    % Background Concentration Gradient: 0 (for no) or 1 (for yes)
electro_phys_Yes = model_Info(17);     % Electrophysiology (FitzHugh-Nagumo): 0 (for no) or 1 (for yes)
d_Springs_Yes = model_Info(18);        % Damped Springs: 0 (for no) or 1 (for yes)
update_D_Springs_Flag = model_Info(19);% Update_Damped_Springs: % 0 (for no) or 1 (for yes)
boussinesq_Yes = model_Info(20);       % Boussinesq Approx.: 0 (for no) or 1 (for yes)
exp_Coeff = model_Info(21);            % Expansion Coefficient (e.g., thermal, etc) for Boussinesq approx.

%Lagrangian Structure Data
ds = Lx / (2*Nx);                   %Lagrangian Spacing
grid_Info(9) = ds;


% Create EULERIAN Mesh (these assume periodicity in x and y)
x = (0:dx:Lx-dx); 
y = (0:dy:Ly-dy)';
%Create x-Mesh
%for i=1:Nx
%    X = [X; x]; 
%end
%Create y-Mesh
%for i=1:Ny
%    Y = [Y y];
%end
[X,Y] = meshgrid(0:dx:Lx-dx,0:dy:Ly-dy);




% % % % % HOPEFULLY WHERE I CAN READ IN INFO!!! % % % % %





% READ IN LAGRANGIAN POINTS %
[Nb,xLag,yLag] = read_Vertex_Points(struct_name);
grid_Info(8) = Nb;          % # Total Number of Lagrangian Pts.
xLag_P = xLag;              % Initialize previous Lagrangian x-Values (for use in muscle-model)
yLag_P = yLag;              % Initialize previous Lagrangian y-Values (for use in muscle-model)

fprintf('\n--> FIBER MODEL INCLUDES: \n');



% READ IN SPRINGS (IF THERE ARE SPRINGS) %
if ( springs_Yes == 1 )
    fprintf('  -Springs and ');
    if update_Springs_Flag == 0
        fprintf('NOT dynamically updating spring properties\n');
    else
        fprintf('dynamically updating spring properties\n');
    end
    springs_info = read_Spring_Points(struct_name);
        %springs_info: col 1: starting spring pt (by lag. discretization)
        %              col 2: ending spring pt. (by lag. discretization)
        %              col 3: spring stiffness
        %              col 4: spring resting lengths
        %              col 5: spring linearity (1=linear, >1 non-linear)
else
    springs_info = 0;  %just to pass placeholder into "please_Find_Lagrangian_Forces_On_Eulerian_grid function"
end






% READ IN BEAMS (IF THERE ARE BEAMS) %
if ( beams_Yes == 1)
    fprintf('  -Beams ("Torsional Springs") and ');
    if update_Beams_Flag == 0
        fprintf('NOT dynamically updating beam properties\n');
    else
        fprintf('dynamically updating beam properties\n');
    end
    
    beams_info = read_Beam_Points(struct_name);
    %beams:      col 1: 1ST PT.
    %            col 2: MIDDLE PT. (where force is exerted)
    %            col 3: 3RD PT.
    %            col 4: beam stiffness
    %            col 5: curavture
else
    beams_info = 0;
end



% READ IN TARGET POINTS (IF THERE ARE TARGET PTS) %
if ( target_pts_Yes == 1)
    fprintf('  -Target Pts. and ');
    if update_Target_Pts == 0
        fprintf('NOT dynamically updating target point properties\n');
    else
        fprintf('dynamically updating target point properties\n');
    end
    
    target_aux = read_Target_Points(struct_name);
    %target_aux: col 1: Lag Pt. ID w/ Associated Target Pt.
    %            col 2: target STIFFNESSES
    target_info(:,1) = target_aux(:,1); %Stores Lag-Pt IDs in col vector
    for i=1:length(target_info(:,1))
        id = target_info(i,1);
        target_info(i,2) = xLag(id);    %Stores Original x-Lags as x-Target Pt. Identities
        target_info(i,3) = yLag(id);    %Stores Original y-Lags as y-Target Pt. Identities
    end
   
    target_info(:,4) = target_aux(:,2); %Stores Target Stiffnesses 
else
    target_info = 0;
end




% READ IN MASS POINTS (IF THERE ARE MASS PTS) %
if ( mass_Yes == 1)
    fprintf('  -Mass Pts. with ');
    if gravity_Yes == 0
        fprintf('NO artificial gravity\n');
    else
        fprintf('artificial gravity\n');
    end
    mass_aux = read_Mass_Points(struct_name);
    %target_aux: col 1: Lag Pt. ID w/ Associated Mass Pt.
    %            col 2: "Mass-spring" stiffness parameter
    %            col 3: "MASS" value parameter
    mass_info(:,1) = mass_aux(:,1); %Stores Lag-Pt IDs in col vector
    for i=1:length(mass_info(:,1))
        id = mass_info(i,1);
        mass_info(i,2) = xLag(id);    %Stores Original x-Lags as x-Mass Pt. Identities
        mass_info(i,3) = yLag(id);    %Stores Original y-Lags as y-Mass Pt. Identities
    end
   
    mass_info(:,4) = mass_aux(:,2);   %Stores "mass-spring" parameter 
    mass_info(:,5) = mass_aux(:,3);   %Stores "MASS" value parameter
else
    mass_info = 0;
end


% CONSTRUCT GRAVITY INFORMATION (IF THERE IS GRAVITY) %
if gravity_Yes == 1
    %gravity_Vec(1) = model_Info(12);     % x-Component of Gravity Vector
    %gravity_Vec(2) = model_Info(13);     % y-Component of Gravity Vector
    xG = model_Info(13);
    yG = model_Info(14);
    normG = sqrt( xG^2 + yG^2 );
    gravity_Info = [gravity_Yes xG/normG yG/normG];
    %   col 1: flag if considering gravity
    %   col 2: x-component of gravity vector (normalized)
    %   col 3: y-component of gravity vector (normalized)
    
    clear xG yG normG;
    
else
    gravity_Info = 0;
end


% READ IN POROUS MEDIA INFO (IF THERE IS POROSITY) %
if ( porous_Yes == 1)
    fprintf('  -Porous Points\n');
    porous_aux = read_Porous_Points(struct_name);
    %porous_aux: col 1: Lag Pt. ID w/ Associated Porous Pt.
    %            col 2: Porosity coefficient
    %            col 3: Flag for derivative stencil!
    porous_info(:,1) = porous_aux(:,1); %Stores Lag-Pt IDs in col vector
    for i=1:length(porous_info(:,1))
        id = porous_info(i,1);
        porous_info(i,2) = xLag(id);    %Stores Original x-Lags as x-Porous Pt. Identities
        porous_info(i,3) = yLag(id);    %Stores Original y-Lags as y-Porous Pt. Identities
    end
   
    porous_info(:,4) = porous_aux(:,2); %Stores Porosity Coefficient 
    porous_info(:,5) = porous_aux(:,3); %Stores flag for derivative stencil
else
    porous_info = 0;
end



% READ IN SPRINGS (IF THERE ARE SPRINGS) %
if ( d_Springs_Yes == 1 )
    fprintf('  -Damped Springs and ');
    if update_D_Springs_Flag == 0
        fprintf('NOT dynamically updating damped spring properties\n');
    else
        fprintf('dynamically updating spring properties\n');
    end
    d_springs_info = read_Damped_Spring_Points(struct_name);
        %springs_info: col 1: starting spring pt (by lag. discretization)
        %              col 2: ending spring pt. (by lag. discretization)
        %              col 3: spring stiffness
        %              col 4: spring resting lengths
        %              col 5: damping coefficients
else
    d_springs_info = 0;  %just to pass placeholder into "please_Find_Lagrangian_Forces_On_Eulerian_grid function"
end








% READ IN MUSCLES (IF THERE ARE MUSCLES) %
if ( muscles_Yes == 1 )
    fprintf('  -MUSCLE MODEL (Force-Velocity / Length-Tension Model)\n');
    muscles_info = read_Muscle_Points(struct_name);
        %         muscles: col 1: MASTER NODE (by lag. discretization)
        %         col 2: SLAVE NODE (by lag. discretization)
        %         col 3: length for max. muscle tension
        %         col 4: muscle constant
        %         col 5: hill parameter, a
        %         col 6: hill parameters, b
        %         col 7: force maximum!
else
    muscles_info = 0;  %just to pass placeholder into "please_Find_Lagrangian_Forces_On_Eulerian_grid function"
end












% READ IN MUSCLES (IF THERE ARE MUSCLES) %
if ( hill_3_muscles_Yes == 1 )
    fprintf('  -MUSCLE MODEL (3 Element Hill Model)\n');
    muscles3_info = read_Hill_3Muscle_Points(struct_name);
        %         muscles: col 1: MASTER NODE (by lag. discretization)
        %         col 2: SLAVE NODE (by lag. discretization)
        %         col 3: length for max. muscle tension
        %         col 4: muscle constant
        %         col 5: hill parameter, a
        %         col 6: hill parameters, b
        %         col 7: force maximum!
        %         col 8: kStiffness of NL Spring
        %         col 9: alpha, degree of non-linearity
else
    muscles3_info = 0;  %just to pass placeholder into "please_Find_Lagrangian_Forces_On_Eulerian_grid function"
end




% SOLVE ELECTROPHYSIOLOGY MODEL FOR PUMPING %
if electro_phys_Yes == 1
    fprintf('  -Electrophysiology model via FitzHugh-Nagumo\n');
    fprintf('\n\n--> Solving Electrophysiology Model...\n');
    [electro_potential, ePhys_Start, ePhys_End] = FitzHugh_Nagumo_1d(dt,T_FINAL);
    ePhys_Ct = 1;
    fprintf('--> Finished Computing Electrophysiology...time for IBM!\n');
else
    electro_potential = 0;
end



fprintf('\n\n--> Background Flow Items\n');
if ( tracers_Yes == 0 ) && (concentration_Yes == 0)
    fprintf('      (No tracers nor other passive scalars immersed in fluid)\n\n');
end

% READ IN TRACERS (IF THERE ARE TRACERS) %
if (tracers_Yes == 1)
   fprintf('  -Tracer Particles included\n');
   [~,xT,yT] = read_Tracer_Points(struct_name);
   tracers = zeros(length(xT),4);
   tracers(1,1) = 1;
   tracers(:,2) = xT;
   tracers(:,3) = yT;
        %tracers_info: col 1: xPt of Tracers
        %              col 2: yPt of Tracers
else
   tracers = 0; 
end


% READ IN CONCENTRATION (IF THERE IS A BACKGROUND CONCENTRATION) %
if ( concentration_Yes == 1 )
    fprintf('  -Background concentration included\n');
    [C,kDiffusion] = read_In_Concentration_Info(struct_name,Nx);
        %C:           Initial background concentration
        %kDiffusion:  Diffusion constant for Advection-Diffusion
else
    C = 0; % placeholder for plotting 
end





% CONSTRUCT BOUSSINESQ INFORMATION (IF USING BOUSSINESQ) %
if boussinesq_Yes == 1
    
    if length(gravity_Info) == 1
        fprintf('\n\n\n READ THE ERROR MESSAGE -> YOU MUST FLAG GRAVITY IN INPUT FILE FOR BOUSSINESQ! :) \n\n\n');
        error('YOU MUST FLAG GRAVITY IN INPUT FILE FOR BOUSSINESQ! :)');
    elseif concentration_Yes == 0
        fprintf('\n\n\n READ THE ERROR MESSAGE -> YOU MUST HAVE BACKGROUND CONCENTRATION FOR BOUSSINESQ! :) \n\n\n');
        error('YOU MUST FLAG CONCENTRATION IN INPUT FILE FOR BOUSSINESQ! :)');
    end
    
    % Forms Boussinesq forcing terms, e.g., (exp_Coeff)*gVector for Momentum Eq.
    [fBouss_X,fBouss_Y] = please_Form_Boussinesq_Forcing_Terms(exp_Coeff,Nx,Ny,gravity_Info);
    
    % Finds initial concentration Laplacian
    Cxx = DD(C,dx,'x');
    Cyy = DD(C,dy,'y');
    laplacian_C = Cxx+Cyy;
    C_0 = zeros(size(C)); % Define background concentration
end







% Initialize the initial velocities to zero.
U = zeros(Ny,Nx);                                % x-Eulerian grid velocity
V = U;                                           % y-Eulerian grid velocity
mVelocity = zeros( length(mass_info(:,1)), 2 );  % mass-Pt velocity 

if arb_ext_force_Yes == 1 
    firstExtForce = 1;                           % initialize external forcing
    indsExtForce = 0;                            % initialize for external forcing computation
end

% ACTUAL TIME-STEPPING IBM SCHEME! 
%(flags for storing structure connects for printing and printing to .vtk)
cter = 0; ctsave = 0; firstPrint = 1; loc = 1; diffy = 1;


% CREATE VIZ_IB2D FOLDER and VISIT FILES
mkdir('viz_IB2d');
mkdir('hier_IB2d_data'); 
%cd('viz_IB2d');
%vizID = fopen('dumps.visit','w'); 
%fprintf(vizID,'!NBLOCKS 6\n');
%cd ..


%Initialize Vorticity, uMagnitude, and Pressure for initial colormap
%Print initializations to .vtk
vort=zeros(Nx,Ny); uMag=vort; p = vort;  lagPts = [xLag yLag zeros(length(xLag),1)]; 
[connectsMat,spacing] = give_Me_Lag_Pt_Connects(ds,xLag,yLag,Nx,springs_Yes,springs_info);
Fxh = vort; Fyh = vort; F_Lag = zeros( length(xLag), 2); 
print_vtk_files(ctsave,vort,uMag,p,U,V,Lx,Ly,Nx,Ny,lagPts,springs_Yes,connectsMat,tracers,concentration_Yes,C,Fxh,Fyh,F_Lag);
fprintf('\n |****** Begin IMMERSED BOUNDARY SIMULATION! ******| \n\n');
fprintf('Current Time(s): %6.6f\n\n',current_time);
ctsave = ctsave+1;


%
%
%
% **************************************************************************
% * * * * * * * * * * * * * BEGIN TIME-STEPPING! * * * * * * * * * * * * * *
% **************************************************************************
%
%
%
while current_time < T_FINAL
   
    
    %
    %
    %**************** Step 1: Update Position of Boundary of membrane at half time-step *******************
    %                           (Variables end with h if it is a half-step)
    %
    [xLag_h, yLag_h] = please_Move_Lagrangian_Point_Positions(U, V, xLag, yLag, xLag, yLag, x, y, dt/2, grid_Info,0);
    
    if mass_Yes == 1
       [mass_info, massLagsOld] = please_Move_Massive_Boundary(dt/2,mass_info,mVelocity); 
    end
    
    if ( ( electro_phys_Yes == 1) && (muscles_Yes == 0) )
        springs_info(ePhys_Start:ePhys_End,3) = ( 1e1*electro_potential(ePhys_Ct,:)') .^4;%( 1e4*electro_potential(ePhys_Ct,:)'.*springs_info(ePhys_Start:ePhys_End,3) ).^4;
        ePhys_Ct = ePhys_Ct + 1;
    end
    
    if ( ( update_Springs_Flag == 1 ) && ( springs_Yes == 1 ) )
       springs_info = update_Springs(dt,current_time,xLag,yLag,springs_info); 
    end
    
    if ( ( update_Target_Pts == 1 ) && ( target_pts_Yes == 1) )
       target_info = update_Target_Point_Positions(dt,current_time,target_info); 
    end
    
    if ( ( update_Beams_Flag == 1 ) && ( beams_Yes == 1) )
       beams_info = update_Beams(dt,current_time,beams_info); 
    end
    
    if ( ( update_D_Springs_Flag == 1 ) && ( d_Springs_Yes == 1) )
       d_springs_info = update_Damped_Springs(dt,current_time,d_springs_info); 
    end
    
    
    %
    %
    %**************** STEP 2: Calculate Force coming from membrane at half time-step ****************
    %
    %
    [Fxh, Fyh, F_Mass_Bnd, F_Lag] =    please_Find_Lagrangian_Forces_On_Eulerian_grid(dt, current_time, xLag_h, yLag_h, xLag_P, yLag_P, x, y, grid_Info, model_Info, springs_info, target_info, beams_info, muscles_info, muscles3_info, mass_info, electro_potential, d_springs_info);
    
    % Once force is calculated, can finish time-step for massive boundary
    if mass_Yes == 1    
        % Update Massive Boundary Velocity
        mVelocity_h = please_Update_Massive_Boundary_Velocity(dt/2,mass_info,mVelocity,F_Mass_Bnd,gravity_Info);
        
        % Update Massive Boundary Position for Time-step
        mass_info(:,[2 3]) = massLagsOld;
        [mass_info,~] = please_Move_Massive_Boundary(dt,mass_info,mVelocity_h); 

        % Update Massive Boundary Velocity for Time-step
        mVelocity = please_Update_Massive_Boundary_Velocity(dt,mass_info,mVelocity,F_Mass_Bnd,gravity_Info); 
    end
    
    % Add artificial force from fluid boundary, if desired. 
    if arb_ext_force_Yes == 1 
        [Fx_Arb, Fy_Arb, firstExtForce, indsExtForce] = please_Compute_External_Forcing(dt, current_time, x, y, grid_Info, U, V, firstExtForce, indsExtForce);
        Fxh = Fxh + Fx_Arb;
        Fyh = Fyh + Fy_Arb;
    end
    
    %
    %
    %***************** STEP 3: Solve for Fluid motion ******************************************
    %
    %
    % Add in effect from BOUSSINESQ
    if boussinesq_Yes == 1
        Fxh = Fxh + rho*fBouss_X*(C);
        Fyh = Fyh + rho*fBouss_Y*(C);
        [Uh, Vh, U, V, p] =   please_Update_Fluid_Velocity(U, V, Fxh, Fyh, rho, mu, grid_Info, dt);
    else
        [Uh, Vh, U, V, p] =   please_Update_Fluid_Velocity(U, V, Fxh, Fyh, rho, mu, grid_Info, dt);
    end
    
    %
    %
    %***************** STEP 4: Update Position of Boundary of membrane again for a half time-step *******
    %
    %
    
    xLag_P = xLag_h;   % Stores old Lagrangian x-Values (for muscle model)
    yLag_P = yLag_h;   % Stores old Lagrangian y-Values (for muscle model)
    %Uh, Vh instead of U,V?
    [xLag, yLag] =     please_Move_Lagrangian_Point_Positions(Uh, Vh, xLag, yLag, xLag_h, yLag_h, x, y, dt, grid_Info,porous_Yes);

    
    %NOTE: SET UP FOR BOTH CLOSED + OPEN SYSTEMS NOW!!!
    if porous_Yes == 1
       [Por_Mat,nX,nY] = please_Compute_Porous_Slip_Velocity(ds,xLag,yLag,porous_info,F_Lag);
       xLag( porous_info(:,1) ) = xLag( porous_info(:,1) ) - dt*( Por_Mat(:,1)+Por_Mat(:,2) ).*nX;
       yLag( porous_info(:,1) ) = yLag( porous_info(:,1) ) - dt*( Por_Mat(:,1)+Por_Mat(:,2) ).*nY;
       porous_info(:,2) = xLag( porous_info(:,1) );    %Stores x-Lags as x-Porous Pt. Identities
       porous_info(:,3) = yLag( porous_info(:,1) );    %Stores y-Lags as y-Porous Pt. Identities
       xLag = mod(xLag, Lx); % If structure goes outside domain
       yLag = mod(yLag, Ly); % If structure goes outside domain
    end
    
    
    % If there are tracers, update tracer positions %
    if tracers_Yes == 1
        %Uh, Vh instead of U,V?
        [xT, yT] = please_Move_Lagrangian_Point_Positions(Uh, Vh, xT, yT, xT, yT, x, y, dt, grid_Info,0); %0 for always no porous tracers
        tracers(:,2) = xT;
        tracers(:,3) = yT;
    end
    
    
    % If there is a background concentration, update concentration-gradient %
    if concentration_Yes == 1
       [C,~] = please_Update_Adv_Diff_Concentration(C,dt,dx,dy,U,V,kDiffusion); 
    end
    
        
    % Plot Lagrangian/Eulerian Dynamics! %
    if ( ( mod(cter,pDump) == 0  ) && ( cter >= pDump ) )
        
        %Compute vorticity, uMagnitude
        vort = give_Me_Vorticity(U,V,dx,dy);
        uMag = give_Me_Magnitude_Velocity(U,V);
        
        %Plot in Matlab
        if pMatlab == 1
            [loc, diffy] = please_Plot_Results(ds,X,Y,U,V,vort,uMag,p,xLag,yLag,lagPlot,velPlot,vortPlot,pressPlot,uMagPlot,firstPrint,loc,diffy,spacing);
        end
        
        %Print .vtk files!
        lagPts = [xLag yLag zeros(length(xLag),1)];
        print_vtk_files(ctsave,vort,uMag',p',U',V',Lx,Ly,Nx,Ny,lagPts,springs_Yes,connectsMat,tracers,concentration_Yes,C,Fxh',Fyh',F_Lag);
        
        %Print Current Time
        fprintf('Current Time(s): %6.6f\n',current_time);
        
        % Prints CFL For Advection:
        maxy = max( max(max(abs(U))), max(max(abs(V))) );
        CFL_adv = (dt/dx)*maxy;
        %CFL_para = dt/dx^2*k;
        %fprintf('CFL: %d\n\n',max(CFL_adv,CFL_para));
        fprintf('CFL: %d\n\n',CFL_adv);
        
        
        %Update print counter for filename index
        ctsave=ctsave+1; firstPrint = 0;
        %pause();
         %electro_potential(ePhys_Ct,:)'
         %springs_info(ePhys_Start:ePhys_End,3)
         %springs_info(:,3)
         %pause();
        
    end

    
    % Update current_time value & counter
    current_time = current_time+dt;
    cter = cter + 1;
    %pause(0.);
    
    
end %ENDS TIME-STEPPING LOOP










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives appropriate string number for filename in printing the
% .vtk files. for viz_IB2d and hier_IB2d_data folders
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_vtk_files(ctsave,vort,uMag,p,U,V,Lx,Ly,Nx,Ny,lagPts,springs_Yes,connectsMat,tracers,concentration_Yes,C,fXGrid,fYGrid,F_Lag)

%Give spacing for grid
dx = Lx/Nx; 
dy = Ly/Ny;

%
% Print all Eulerian Data and Lagrangian Position Data to viz_IB2d folder
%
cd('viz_IB2d'); %Go into viz_IB2d directory

    %Find string number for storing files
    strNUM = give_String_Number_For_VTK(ctsave);
    vortfName = ['Omega.' strNUM '.vtk'];
    uMagfName = ['uMag.' strNUM '.vtk'];
    pfName = ['P.' strNUM '.vtk'];
    uXName = ['uX.' strNUM '.vtk'];
    uYName = ['uY.' strNUM '.vtk'];
    fXName = ['fX.' strNUM '.vtk'];
    fYName = ['fY.' strNUM '.vtk'];
    fMagName = ['fMag.' strNUM '.vtk'];
    velocityName = ['u.' strNUM '.vtk'];
    lagPtsName = ['lagsPts.' strNUM '.vtk'];

    %Print Lagrangian Pts to .vtk format
    savevtk_points(lagPts, lagPtsName, 'lagPts');

    % Print Spring Connections (* if springs *)
    if springs_Yes == 1
         %Print Lagrangian Pts w/ CONNECTIONS to .vtk format
        lagPtsConName=['lagPtsConnect.' strNUM '.vtk'];
        savevtk_points_connects(lagPts, lagPtsConName, 'lagPtsConnected',connectsMat);
    end
    
    %Print Tracer Pts (*if tracers*)
    if tracers(1,1) == 1
        tracersPtsName = ['tracer.' strNUM '.vtk'];
        %tMatrix = tracers(:,2:4);
        savevtk_points(tracers(:,2:4),tracersPtsName, 'tracers'); 
    end

    %Print another cycle to .visit file
    %fprintf(vizID,[vortfName '\n']);
    %fprintf(vizID,[uMagfName '\n']);
    %fprintf(vizID,[pfName '\n']);
    %fprintf(vizID,[uXName '\n']);
    %fprintf(vizID,[uYName '\n']);
    %fprintf(vizID,[velocityName '\n']);


    %Print SCALAR DATA (i.e., colormap data) to .vtk file
    savevtk_scalar(vort, vortfName, 'Omega',dx,dy);
    savevtk_scalar(uMag, uMagfName, 'uMag',dx,dy);
    %savevtk_scalar(p, pfName, 'P',dx,dy);
    %savevtk_scalar(U, uXName, 'uX',dx,dy);
    %savevtk_scalar(V, uYName, 'uY',dx,dy);
    %savevtk_scalar(fXGrid, fXName, 'fX',dx,dy);
    %savevtk_scalar(fYGrid, fYName, 'fY',dx,dy);
    %savevtk_scalar(sqrt( fXGrid.^2 + fYGrid.^2 ), fMagName, 'fMag',dx,dy);


    if concentration_Yes == 1
        confName = ['concentration.' strNUM '.vtk'];
        savevtk_scalar(C', confName, 'Concentration',dx,dy);
    end


    %Print VECTOR DATA (i.e., velocity data) to .vtk file
    savevtk_vector(U, V, velocityName, 'u',dx,dy)

%Get out of viz_IB2d folder
cd ..

if length( lagPts ) <= 5
    %
    % Print Lagrangian Force Data to hier_IB2d_data folder (if <= 5 lag pts)
    %
    cd('hier_IB2d_data'); %change directory to hier-data folder
    
        % Save x-y force data!
        fLag_XName = ['fX_Lag.' strNUM '.vtk'];
        fLag_YName = ['fY_Lag.' strNUM '.vtk'];
        savevtk_points_with_scalar_data( lagPts, F_Lag(:,1), fLag_XName, 'fX_Lag');
        savevtk_points_with_scalar_data( lagPts, F_Lag(:,2), fLag_YName, 'fY_Lag');

        % Save force Magnitude (no normal/tangential -> not enough points)
        fMagName = ['fMag.' strNUM '.vtk'];
        fLagMag = sqrt( F_Lag(:,1).^2 + F_Lag(:,2).^2 ); % Compute magnitude of forces on boundary
        savevtk_points_with_scalar_data( lagPts, fLagMag, fMagName, 'fMag');
    cd ..
else

    %
    % Print Lagrangian Force Data to hier_IB2d_data folder
    %
    [F_Tan_Mag,F_Normal_Mag] = please_Compute_Normal_Tangential_Forces_On_Lag_Pts(lagPts,F_Lag);
    %
    cd('hier_IB2d_data'); %change directory to hier-data folder
        
        % Save x-y force data!
        fLag_XName = ['fX_Lag.' strNUM '.vtk'];
        fLag_YName = ['fY_Lag.' strNUM '.vtk'];
        savevtk_points_with_scalar_data( lagPts, F_Lag(:,1), fLag_XName, 'fX_Lag');
        savevtk_points_with_scalar_data( lagPts, F_Lag(:,2), fLag_YName, 'fY_Lag');
    
        % Save force Magnitude, Mag. Normal, Mag. Tangential
        fMagName = ['fMag.' strNUM '.vtk'];
        fNormalName = ['fNorm.' strNUM '.vtk'];
        fTangentName = ['fTan.' strNUM '.vtk'];

        fLagMag = sqrt( F_Lag(:,1).^2 + F_Lag(:,2).^2 ); % Compute magnitude of forces on boundary

        % savevtk_points_with_scalar_data( lagPts, fLagMag, fMagName, 'fMag');
        %savevtk_points_with_scalar_data( lagPts, F_Normal_Mag, fNormalName, 'fNorm');
        %savevtk_points_with_scalar_data( lagPts, F_Tan_Mag, fTangentName, 'fTan');

    cd .. % Get out of hier_IB2d_data folder

end % ENDS IF-STATEMENT


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: give me Connects Vector for printing Lagrangian .vtk info!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [connectsMat,space] = give_Me_Lag_Pt_Connects(ds,xLag,yLag,Nx,springs_Yes,springs_info)

        %springs_info: col 1: starting spring pt (by lag. discretization)
        %              col 2: ending spring pt. (by lag. discretization)
        %              col 3: spring stiffness
        %              col 4: spring resting lengths
        %              col 5: spring linearity (1=linear, >1 non-linear)



space = 20*ds;

if springs_Yes == 1
    connectsMat(:,1) = springs_info(:,1)-1; % -1 For Cpp notation (and .vtk counting)
    connectsMat(:,2) = springs_info(:,2)-1; % -1 For Cpp notation (and .vtk counting)
else
    
    connectsMat = 0;
end

%N = length(xLag);

% if Nx <= 32
%     space = 5*ds;
% elseif Nx <= 64 
%    space = 5*ds; 
% elseif Nx <=128
%    space = 5*ds;
% elseif Nx <=256
%     space = 10*ds;
% elseif Nx <= 512
%     space = 20*ds;
% else
%     space = 40*ds;
% end
    

% ct = 1;
% for i=1:N
%     if i<N
%         x1=xLag(i); x2=xLag(i+1);
%         y1=yLag(i); y2=yLag(i+1);
%         dist = sqrt( (x1-x2)^2 + (y1-y2)^2 );
%         if dist < space
%             connectsMat(ct,1) = i-1; %For Cpp notation (and .vtk counting)
%             connectsMat(ct,2) = i;   %For Cpp notation (and .vtk counting)
%         ct=ct+1;
%         end
%     elseif i==N
%         x1=xLag(N); x2=xLag(1);
%         y1=yLag(N); y2=yLag(1);
%         dist = sqrt( (x1-x2)^2 + (y1-y2)^2 );
%         if dist < space
%             connectsMat(ct,1) = N-1; %For Cpp notation (and .vtk counting)
%             connectsMat(ct,2) = 0;   %For Cpp notation (and .vtk counting)
%         ct=ct+1;
%         end
%     end
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives appropriate string number for filename in printing the
% .vtk files.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function strNUM = give_String_Number_For_VTK(num)

%num: # of file to be printed

if num < 10
    strNUM = ['000' num2str(num)];
elseif num < 100
    strNUM = ['00' num2str(num)];
elseif num<1000
    strNUM = ['0' num2str(num)];
else
    strNUM = num2str(num);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Computes vorticity from two matrices, U and V, where each
% matrix is the discretized field of velocity values either for x or y,
% respectively.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vort = give_Me_Vorticity(U,V,dx,dy)

% w = ( dv/dx - du/dy )\hat{z}

%Compute dv/dx using central differencing! (maintains periodicity)
dvdx = D(V,dx,'x');

%Compute du/dy using central differencing! (maintains periodicity)
dudy = D(U,dy,'y');

%Compute vorticity
vort = ( dvdx - dudy );

%Take transpose so all lines up
vort = vort';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Computes vorticity from two matrices, U and V, where each
% matrix is the discretized field of velocity values either for x or y,
% respectively.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function uMag = give_Me_Magnitude_Velocity(U,V)

% U: x-directed velocity
% V: y-directed velocity

% Compute magnitude of velocity
uMag = ( U.^2 + V.^2 ).^(1/2);
 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints matrix vector data to vtk formated file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function savevtk_points_connects( X, filename, vectorName,connectsMat)

%X is matrix of size Nx3

N = length( X(:,1) );
Nc = length( connectsMat(:,1) );

%TRY PRINTING THEM AS UNSTRUCTURED_GRID
file = fopen (filename, 'w');
fprintf(file, '# vtk DataFile Version 2.0\n');
fprintf(file, [vectorName '\n']);
fprintf(file, 'ASCII\n');
fprintf(file, 'DATASET UNSTRUCTURED_GRID\n\n');
%
fprintf(file, 'POINTS %i float\n', N);
for i=1:N
    fprintf(file, '%.15e %.15e %.15e\n', X(i,1),X(i,2),X(i,3));
end
fprintf(file,'\n');
%
fprintf(file,'CELLS %i %i\n',Nc,3*Nc); %First: # of "Cells", Second: Total # of info inputed following
for s=1:Nc
    fprintf(file,'%i %i %i\n',2, connectsMat(s,1), connectsMat(s,2) );
end
fprintf(file,'\n');
%
fprintf(file,'CELL_TYPES %i\n',Nc); % N = # of "Cells"
for i=1:Nc
   fprintf(file,'3 '); 
end
fprintf(file,'\n');
fclose(file);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints matrix vector data to vtk formated file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function savevtk_points( X, filename, vectorName)

%X is matrix of size Nx3

N = length( X(:,1) );


%TRY PRINTING THEM AS UNSTRUCTURED_GRID
file = fopen (filename, 'w');
fprintf(file, '# vtk DataFile Version 2.0\n');
fprintf(file, [vectorName '\n']);
fprintf(file, 'ASCII\n');
fprintf(file, 'DATASET UNSTRUCTURED_GRID\n\n');
%
fprintf(file, 'POINTS %i float\n', N);
for i=1:N
    fprintf(file, '%.15e %.15e %.15e\n', X(i,1),X(i,2),X(i,3));
end
fprintf(file,'\n');
%
fprintf(file,'CELLS %i %i\n',N,2*N); %First: # of "Cells", Second: Total # of info inputed following
for s=0:N-1
    fprintf(file,'%i %i\n',1,s);
end
fprintf(file,'\n');
%
fprintf(file,'CELL_TYPES %i\n',N); % N = # of "Cells"
for i=1:N
   fprintf(file,'1 '); 
end
fprintf(file,'\n');
fclose(file);



%TRY PRINTING THEM AS POLYGONAL DATA
% file = fopen (filename, 'w');
% fprintf(file, '# vtk DataFile Version 2.0\n');
% fprintf(file, [vectorName '\n']);
% fprintf(file, 'ASCII\n');
% fprintf(file, 'DATASET STRUCTURED_GRID\n');
% fprintf(file, 'DIMENSIONS 64 1 1\n');
% fprintf(file, 'POINTS %i float\n', N);
% for i=1:N
% fprintf(file, '%.15e %.15e %.15e\n', X(i,1),X(i,2),X(i,3));
% end
% fprintf(file,'1.1 1.1 0\n');
% fprintf(file,'CELL_DATA 1\n');
% fprintf(file,'POINT_DATA %u \n',N);
% fprintf(file,'FIELD FieldData 1\n');
% fprintf(file,'nodal 1 %i float\n',N);
% fprintf(file,'0 1 1.1 2\n');
% fprintf(file,'SCALARS nodal float\n');
% fprintf(file,['SCALARS ' vectorName ' float 1 \n']);
% fprintf(file,'LOOKUP_TABLE default\n');


% TRY PRINTING THEM AS POINTS
% file = fopen (filename, 'w');
% fprintf(file, '# vtk DataFile Version 2.0\n');
% fprintf(file, 'Cube example\n');
% fprintf(file, 'ASCII\n');
% fprintf(file, 'DATASET UNSTRUCTURED_GRID\n');
% fprintf(file, 'POINTS %i float\n', N);
% for i=1:N
% fprintf(file, '%.15e %.15e %.15e\n', X(i,1),X(i,2),X(i,3));
% end
% fprintf(file,'POINT_DATA %u \n',N);
% fprintf(file,['SCALARS ' vectorName ' float 1 \n']);
% fprintf(file,'LOOKUP_TABLE default\n');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints matrix vector data to vtk formated file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function savevtk_vector(X, Y, filename, vectorName,dx,dy)
%  savevtkvector Save a 3-D vector array in VTK format
%  savevtkvector(X,Y,Z,filename) saves a 3-D vector of any size to
%  filename in VTK format. X, Y and Z should be arrays of the same
%  size, each storing speeds in the a single Cartesian directions.
    if (size(X) ~= size(Y))
        fprint('Error: velocity arrays of unequal size\n'); return;
    end
    [nx, ny, nz] = size(X);
    fid = fopen(filename, 'wt');
    fprintf(fid, '# vtk DataFile Version 2.0\n');
    fprintf(fid, 'Comment goes here\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, '\n');
    fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
    fprintf(fid, 'DIMENSIONS    %d   %d   %d\n', nx, ny, nz);
    fprintf(fid, '\n');
    fprintf(fid, 'ORIGIN    0.000   0.000   0.000\n');
    %fprintf(fid, 'SPACING   1.000   1.000   1.000\n'); if want [1,32]x[1,32] rather than [0,Lx]x[0,Ly]
    fprintf(fid, ['SPACING   ' num2str(dx) ' '  num2str(dy) '   1.000\n']);
    fprintf(fid, '\n');
    fprintf(fid, 'POINT_DATA   %d\n', nx*ny);
    fprintf(fid, ['VECTORS ' vectorName ' double\n']);
    fprintf(fid, '\n');
    for a=1:nz
        for b=1:ny
            for c=1:nx
                fprintf(fid, '%f ', X(c,b,1));
                fprintf(fid, '%f ', Y(c,b,1));
                fprintf(fid, '%f ', 0);
            end
            fprintf(fid, '\n');
        end
    end
    fclose(fid);
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints scalar matrix to vtk formated file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function savevtk_scalar(array, filename, colorMap,dx,dy)
%  savevtk Save a 3-D scalar array in VTK format.
%  savevtk(array, filename) saves a 3-D array of any size to
%  filename in VTK format.
    [ny, nx, nz] = size(array);
    fid = fopen(filename, 'wt');
    fprintf(fid, '# vtk DataFile Version 2.0\n');
    fprintf(fid, 'Comment goes here\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, '\n');
    fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
    fprintf(fid, 'DIMENSIONS    %d   %d   %d\n', ny, nx, nz);
    fprintf(fid, '\n');
    fprintf(fid, 'ORIGIN    0.000   0.000   0.000\n');
    %fprintf(fid, 'SPACING   1.000   1.000   1.000\n'); if want [1,32]x[1,32] rather than [0,Lx]x[0,Ly]
    fprintf(fid, ['SPACING   ' num2str(dx) ' '   num2str(dy) '   1.000\n']);
    fprintf(fid, '\n');
    fprintf(fid, 'POINT_DATA   %d\n', nx*ny*nz);
    fprintf(fid, ['SCALARS ' colorMap ' double\n']);
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid, '\n');
    for a=1:nz
        for b=1:nx
            for c=1:ny
                fprintf(fid, '%d ', array(c,b,a));
            end
            fprintf(fid, '\n');
        end
    end
    fclose(fid);
return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints Lagrangian pt data w/ associated scalar to vtk formated file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function savevtk_points_with_scalar_data( X, scalarArray, filename, vectorName)

%X is matrix of size Nx3

N = length( X(:,1) );


%TRY PRINTING THEM AS UNSTRUCTURED_GRID
file = fopen (filename, 'w');
fprintf(file, '# vtk DataFile Version 2.0\n');
fprintf(file, [vectorName '\n']);
fprintf(file, 'ASCII\n');
fprintf(file, 'DATASET UNSTRUCTURED_GRID\n\n');
%
fprintf(file, 'POINTS %i float\n', N);
for i=1:N
    fprintf(file, '%.15e %.15e %.15e\n', X(i,1),X(i,2),X(i,3));
end
fprintf(file,'\n');
%
fprintf(file, 'POINT_DATA   %d\n', N);
fprintf(file, ['SCALARS ' vectorName ' double\n']);
fprintf(file, 'LOOKUP_TABLE default\n');
fprintf(file, '\n');
    for i=1:N
        fprintf(file, '%d ', scalarArray(i,1));
        fprintf(file, '\n');
    end

fclose(file);
    
    
%     for a=1:nz
%         for b=1:ny
%             for c=1:nx
%                 fprintf(fid, '%d ', array(c,b,a));
%             end
%             fprintf(fid, '\n');
%         end
%     end
    



% fprintf(file,'CELLS %i %i\n',N,2*N); %First: # of "Cells", Second: Total # of info inputed following
% for s=0:N-1
%     fprintf(file,'%i %i\n',1,s);
% end
% fprintf(file,'\n');
% %
% fprintf(file,'CELL_TYPES %i\n',N); % N = # of "Cells"
% for i=1:N
%    fprintf(file,'1 '); 
% end
% fprintf(file,'\n');
% fclose(file);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in the # of vertex pts and all the vertex pts from the
%           .vertex file.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [N,xLag,yLag] = read_Vertex_Points(struct_name)

filename = [struct_name '.vertex'];  %Name of file to read in
fileID = fopen(filename);

% Read in the file, use 'CollectOutput' to gather all similar data together
% and 'CommentStyle' to to end and be able to skip lines in file.
C = textscan(fileID,'%f %f','CollectOutput',1);


fclose(fileID);     %Close the data file.

vertices = C{1};    %Stores all read in data in vertices (N+1,2) array

N = vertices(1,1);  % # of Lagrangian Pts
xLag = zeros(N,1);  % Initialize storage for Lagrangian Pts.
yLag = xLag;        % Initialize storage for Lagrangian Pts.

for i=1:N
   xLag(i,1) = vertices(i+1,1); %Stores x-values of Lagrangian Mesh
   yLag(i,1) = vertices(i+1,2); %Stores y-values of Lagrangian Mesh
   
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in the # of tracer pts and all the tracer pts from the
%           .tracer file.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [N,xLag,yLag] = read_Tracer_Points(struct_name)

filename = [struct_name '.tracer'];  %Name of file to read in
fileID = fopen(filename);

% Read in the file, use 'CollectOutput' to gather all similar data together
% and 'CommentStyle' to to end and be able to skip lines in file.
C = textscan(fileID,'%f %f','CollectOutput',1);


fclose(fileID);     %Close the data file.

tracers = C{1};    %Stores all read in data in vertices (N+1,2) array

N = tracers(1,1);  % # of Lagrangian Pts
xLag = zeros(N,1);  % Initialize storage for Lagrangian Pts.
yLag = xLag;        % Initialize storage for Lagrangian Pts.

for i=1:N
   xLag(i,1) = tracers(i+1,1); %Stores x-values of Tracer Lagrangian Mesh
   yLag(i,1) = tracers(i+1,2); %Stores y-values of Tracer Lagrangian Mesh
   
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in the diffusion coefficient and initial concentration, C
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C,kDiff] = read_In_Concentration_Info(struct_name,N)

filename = [struct_name '.concentration'];  %Name of file to read in

strstr = '%f';
for i=1:N-1
   strstr = [strstr ' %f'];
end

fileID = fopen(filename);

    % Read in the file, use 'CollectOutput' to gather all similar data together
    % and 'CommentStyle' to to end and be able to skip lines in file.
    C = textscan(fileID,strstr,'CollectOutput',1);

fclose(fileID);        %Close the data file.

con_info = C{1};    %Stores all read in data 

%Store all elements on .concentration file 
kDiff = con_info(1,1);     %coefficient of diffusion
C = con_info(2:end,1:end); %initial concentration







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in the # of springs and all MASTER NODEs, SLAVE NODEs,
%           spring STIFFNESSES, spring RESTING LENGTHS, spring linearities
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function springs = read_Spring_Points(struct_name)

filename = [struct_name '.spring'];  %Name of file to read in

fileID = fopen(filename);

    % Read in the file, use 'CollectOutput' to gather all similar data together
    % and 'CommentStyle' to to end and be able to skip lines in file.
    C = textscan(fileID,'%f %f %f %f %f','CollectOutput',1);
    
fclose(fileID);        %Close the data file.

spring_info = C{1};    %Stores all read in data in vertices (N+1,2) array

%Store all elements on .spring file into a matrix starting w/ 2nd row of read in data to 4th col.
springs = spring_info(2:end,1:4);

%Stores last column to check for non-linearities
NLspr = spring_info(2:end,5);

%Finds indices of non-linearities
inds = isnan(NLspr);
springs2 = zeros(length(inds),1);

%Stores linearity coefficients
for i=1:length(inds)
    if inds(i) == 1
        springs2(i) = 1;
    else
        springs2(i) = NLspr(i);
    end
end

springs = [springs springs2];

clear springs2 inds NPspr spring_info;

%springs: col 1: starting spring pt (by lag. discretization)
%         col 2: ending spring pt. (by lag. discretization)
%         col 3: spring stiffness
%         col 4: spring resting lengths
%         col 5: linearity of spring (1=linear, > 1 non-linear spring)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in the # of TARGET PTS, TARGET-PT-NODEs, and their
%           Target-STIFFNESSES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function masses = read_Mass_Points(struct_name)

filename = [struct_name '.mass'];  %Name of file to read in

fileID = fopen(filename);

    % Read in the file, use 'CollectOutput' to gather all similar data together
    % and 'CommentStyle' to to end and be able to skip lines in file.
    C = textscan(fileID,'%f %f %f','CollectOutput',1);

fclose(fileID);      %Close the data file.

mass_info = C{1};    %Stores all read in data in vertices (N+1,2) array

%Store all elements on .mass file into a matrix starting w/ 2nd row of read in data.
masses = mass_info(2:end,1:3);

%masses:  col 1: Lag Pt. ID w/ Associated Mass Pt.
%         col 2: "Mass-Spring" stiffness Parameter
%         col 3: Mass Value Parameter





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in the # of TARGET PTS, TARGET-PT-NODEs, and their
%           Target-STIFFNESSES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function targets = read_Target_Points(struct_name)

filename = [struct_name '.target'];  %Name of file to read in

fileID = fopen(filename);

    % Read in the file, use 'CollectOutput' to gather all similar data together
    % and 'CommentStyle' to to end and be able to skip lines in file.
    C = textscan(fileID,'%f %f','CollectOutput',1);

fclose(fileID);        %Close the data file.

targets_info = C{1};    %Stores all read in data in vertices (N+1,2) array

%Store all elements on .spring file into a matrix starting w/ 2nd row of read in data.
targets = targets_info(2:end,1:2);

%targets: col 1: Lag Pt. ID w/ Associated Target Pt.
%         col 2: target STIFFNESSES







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in the # of POROUS PTS, POROUS-PT-NODEs, and their
%           POROUSITY-COEFFICIENTS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function porosity = read_Porous_Points(struct_name)

filename = [struct_name '.porous'];  %Name of file to read in

fileID = fopen(filename);

    % Read in the file, use 'CollectOutput' to gather all similar data together
    % and 'CommentStyle' to to end and be able to skip lines in file.
    C = textscan(fileID,'%f %f %f','CollectOutput',1);

fclose(fileID);        %Close the data file.

porous_info = C{1};    %Stores all read in data in vertices (N+1,2) array

%Store all elements on .spring file into a matrix starting w/ 2nd row of read in data.
porosity = porous_info(2:end,1:3);

%porous:  col 1: Lag Pt. ID w/ Associated Porous Pt.
%         col 2: Porosity coefficient





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in the # of beams and all 1st Pt, MIDDLE Pt, and 3rd Pt
%           beam STIFFNESSES, and CURVATURE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function beams = read_Beam_Points(struct_name)

filename = [struct_name '.beam'];  %Name of file to read in

fileID = fopen(filename);

    % Read in the file, use 'CollectOutput' to gather all similar data together
    % and 'CommentStyle' to to end and be able to skip lines in file.
    C = textscan(fileID,'%f %f %f %f %f','CollectOutput',1);

fclose(fileID);      %Close the data file.

beam_info = C{1};    %Stores all read in data in vertices (N+1,2) array

%Store all elements on .beam file into a matrix starting w/ 2nd row of read in data.
beams = beam_info(2:end,1:5);

    %beams:      col 1: 1ST PT.
    %            col 2: MIDDLE PT. (where force is exerted)
    %            col 3: 3RD PT.
    %            col 4: beam stiffness
    %            col 5: curavture



    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in the # of beams and all 1st Pt, MIDDLE Pt, and 3rd Pt
%           beam STIFFNESSES, and CURVATURE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dSprings = read_Damped_Spring_Points(struct_name)

filename = [struct_name '.d_spring'];  %Name of file to read in

fileID = fopen(filename);

    % Read in the file, use 'CollectOutput' to gather all similar data together
    % and 'CommentStyle' to to end and be able to skip lines in file.
    C = textscan(fileID,'%f %f %f %f %f','CollectOutput',1);

fclose(fileID);      %Close the data file.

dSprings_info = C{1};    %Stores all read in data in vertices (N+1,2) array

%Store all elements on .d_springs file into a matrix starting w/ 2nd row of read in data.
dSprings = dSprings_info(2:end,1:5);

% dSprings: col 1: starting spring pt (by lag. discretization)
%           col 2: ending spring pt. (by lag. discretization)
%           col 3: spring stiffness
%           col 4: spring resting lengths
%           col 5: damping coefficient
    
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in the # of muscles and all MASTER NODEs, SLAVE NODEs,
%           length for max. muscle tension, muscle constant, hill
%           parameters (a and b), and Force-Max
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function muscles = read_Muscle_Points(struct_name)

filename = [struct_name '.muscle'];  %Name of file to read in

fileID = fopen(filename);

    % Read in the file, use 'CollectOutput' to gather all similar data together
    % and 'CommentStyle' to to end and be able to skip lines in file.
    C = textscan(fileID,'%f %f %f %f %f %f %f','CollectOutput',1);

fclose(fileID);        %Close the data file.

muscle_info = C{1};    %Stores all read in data in vertices (N+1,2) array

%Store all elements on .muscle file into a matrix starting w/ 2nd row of read in data.
muscles = muscle_info(2:end,1:7);

%muscles: col 1: MASTER NODE (by lag. discretization)
%         col 2: SLAVE NODE (by lag. discretization)
%         col 3: length for max. muscle tension
%         col 4: muscle constant
%         col 5: hill parameter, a
%         col 6: hill parameters, b
%         col 7: force maximum!
    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in the # of muscles and all MASTER NODEs, SLAVE NODEs,
%           length for max. muscle tension, muscle constant, hill
%           parameters (a and b), and Force-Max
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function muscles = read_Hill_3Muscle_Points(struct_name)

filename = [struct_name '.3_muscle'];  %Name of file to read in

fileID = fopen(filename);

    % Read in the file, use 'CollectOutput' to gather all similar data together
    % and 'CommentStyle' to to end and be able to skip lines in file.
    C = textscan(fileID,'%f %f %f %f %f %f %f %f %f','CollectOutput',1);

fclose(fileID);        %Close the data file.

muscle_info = C{1};    %Stores all read in data in vertices (N+1,2) array

%Store all elements on .muscle file into a matrix starting w/ 2nd row of read in data.
muscles = muscle_info(2:end,1:9);

%muscles: col 1: MASTER NODE (by lag. discretization)
%         col 2: SLAVE NODE (by lag. discretization)
%         col 3: length for max. muscle tension
%         col 4: muscle constant
%         col 5: hill parameter, a
%         col 6: hill parameters, b
%         col 7: force maximum!
%         col 8: NL Spring stiffness, kSpr
%         col 9: NL Spring deg. of non-linearity, alpha


    