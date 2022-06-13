%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  nickabattista@gmail.com
% Date Created: May 27th, 2015
% Institution Created: UNC-CH
% Current Institution: The College of New Jersey (TCNJ)
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
% If you would like us to add a specific muscle model, please let Nick (nickabattista@gmail.com) know.
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Actual DRIVER of the code, where the time-stepping occurs ->
%           gets called by main2d to do the "magic" :)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,Y,U,V,xLag,yLag,C] = IBM_Driver(Fluid_Params,Grid_Params,Time_Params,Lag_Struct_Params,Output_Params,Lag_Name_Params,Con_Params)


%
%    2D IMMERSED BOUNDARY SOLVER ON RECTANGULAR DOMAIN w/ PERIODIC BOUNDARIES
%    
%    x-Momentum Conservation: rho*u_t = -rho*u*u_x + rho*v*u_y + mu*laplacian(u) - p_x + F_x
%    y-Momentum 
%
%    Incompressibility: u_x + v_y = 0
%
%    LagPts/dt = int{ u(x,t) delta( x - LagPts(s,t) ) dx }
%    F_x = int{ fx(s,t) delta(x - LagPts(s,t)) ds }
%    F_y = int{ fy(s,t) delta(x - LagPts(s,t)) ds }
%
fprintf('\n________________________________________________________________________________\n\n');
fprintf('\n---------------->>                 IB2d                      <<----------------\n');
fprintf('\n________________________________________________________________________________\n\n');
fprintf('If using the code for research purposes please cite the following three papers: \n');
fprintf('     [1] N.A. Battista, A.J. Baird, L.A. Miller, A mathematical model and MATLAB code for muscle-fluid-structure simulations, Integ. Comp. Biol. 55(5):901-911 (2015)\n');
fprintf('     [2] N.A. Battista, W.C. Strickland, L.A. Miller, IB2d a Python and MATLAB implementation of the immersed boundary method, Bioinspir. Biomim. 12(3):036003 (2017)\n');
fprintf('     [3] N.A. Battista, W.C. Strickland, A. Barrett, L.A. Miller, IB2d Reloaded: A more powerful Python and MATLAB implementation of the immersed boundary method, Math. Method. Appl. Sci. 41(18):8455-8480 (2018).');
fprintf('\n________________________________________________________________________________');

fprintf('\n\n\n |****** Prepping Immersed Boundary Simulation ******|\n');
fprintf('\n\n--> Reading input data for simulation...\n\n');

%
% ** IBM_DRIVER INPUT DEFINITIONS ** :
%
%               Fluid_Params(1): mu
%                           (2): density
%
%               Grid_Params(1): Nx
%                          (2): Ny
%                          (3): Lx
%                          (4): Ly
%                          (5): Supp
%
%               Time_Params(1): Tfinal (end time of simulation)
%                          (2): dt (time-step)
%
%               Lag_Struct_Params(1): springs
%                                (2): update_springs
%                                (3): target points
%                                (4): update_target_points
%                                (5): beams (torsional beams)
%                                (6): update_beams
%                                 .         .
%                                 .         .
%                                 .         .
%
%               Output_Params(1):  print_dump
%                            (2):  plot_Matlab
%                            (3):  plot_LagPts
%                            (4):  plot_Velocity
%                            (5):  plot_Vorticity
%                            (6):  plot_MagVelocity
%                            (7):  plot_Pressure
%                            (8):  plot_Concentration
%                            (9):  save_Vorticity 
%                            (10):  save_Pressure 
%                            (11): save_uVec 
%                            (12): save_uMag 
%                            (13): save_uX 
%                            (14): save_uY 
%                            (15): save_fMag 
%                            (16): save_fX 
%                            (17): save_fY 
%                            (18): save_hier 
%
%               Con_Params(1): concentration
%                         (2): kDiff
%                         (3): advection
%                         (4): source
%                         (5): k_source
%                         (6): c_inf
%                         (7): num_con

% SIMULATION NAME STRING TO RUN .vertex, .spring, etc. %
struct_name = char(Lag_Name_Params);

% FLUID PARAMETER VALUES STORED %
mu = Fluid_Params(1);      % Dynamic Viscosity
rho = Fluid_Params(2);     % Density

% TEMPORAL INFORMATION VALUES STORED %
T_FINAL = Time_Params(1);      % Final simulation time
dt = Time_Params(2);           % Time-step
restart_Flag = Time_Params(3); % Restart Flag (1=YES, RESTART FROM PREVIOUS DATA, 0=No)
current_time = 0.0;            % initialize start of simulation to time, 0 

% GRID INFO %
Nx =   Grid_Params(1);                % # of Eulerian pts. in x-direction
Ny =   Grid_Params(2);                % # of Eulerian pts. in y-direction
Lx =   Grid_Params(3);                % Length of Eulerian grid in x-coordinate
Ly =   Grid_Params(4);                % Length of Eulerian grid in y-coordinate
dx =   Grid_Params(3)/Grid_Params(1); % Spatial-size in x
dy =   Grid_Params(4)/Grid_Params(2); % Spatial-size in y
supp = Grid_Params(5);                % Delta-function support
grid_Info = [Nx Ny Lx Ly dx dy supp]; % Store for passing values into IB time-loop sub-functions
                                      % NOTE: grid_Info(8) = Nb, grid_Info(9) = ds [THEY GET STORED LATER]

% PRINTING/PLOTTING INFO %
pDump = Output_Params(1);          % Print (Plot) Dump interval
pMatlab = Output_Params(2);        % Plot in Matlab? (1=YES,0=NO)
lagPlot = Output_Params(3);        % Plot LAGRANGIAN PTs ONLY in Matlab
velPlot = Output_Params(4);        % Plot LAGRANGIAN PTs + VELOCITY FIELD in Matlab
vortPlot = Output_Params(5);       % Plot LAGRANGIAN PTs + VORTICITY colormap in Matlab
uMagPlot = Output_Params(6);       % Plot LAGRANGIAN PTs + MAGNITUDE OF VELOCITY colormap in Matlab
pressPlot = Output_Params(7);      % Plot LAGRANGIAN PTs + PRESSURE colormap in Matlab
conPlot = Output_Params(8);      % Plot LAGRANGIAN PTs + PRESSURE colormap in Matlab


% MODEL STRUCTURE DATA STORED %
springs_Yes = Lag_Struct_Params(1);           % Springs: 0 (for no) or 1 (for yes) 
update_Springs_Flag = Lag_Struct_Params(2);   % Update_Springs: 0 (for no) or 1 (for yes)
target_pts_Yes = Lag_Struct_Params(3);        % Target_Pts: 0 (for no) or 1 (for yes)
update_Target_Pts = Lag_Struct_Params(4);     % Update_Target_Pts: 0 (for no) or 1 (for yes)
beams_Yes = Lag_Struct_Params(5);             % Beams: 0 (for no) or 1 (for yes)
update_Beams_Flag = Lag_Struct_Params(6);     % Update_Beams: 0 (for no) or 1 (for yes)
nonInv_beams_Yes = Lag_Struct_Params(7);      % Beams (non-invariant): 0 (for no) or 1 (for yes)
update_nonInv_Beams_Flag = Lag_Struct_Params(8); % Update_nonInv_Beams: 0 (for no) or 1 (for yes)
muscles_Yes = Lag_Struct_Params(9);           % FV-LT Muscles: 0 (for no) or 1 (for yes)
hill_3_muscles_Yes = Lag_Struct_Params(10);    % Hill 3-Element Muscle: 0 (for no) or 1 (for yes)
arb_ext_force_Yes = Lag_Struct_Params(11);     % Arbitrary External Force: 0 (for no) or 1 (for yes)
tracers_Yes = Lag_Struct_Params(12);          % Tracers: 0 (for no) or 1 (for yes)
mass_Yes = Lag_Struct_Params(13);             % Mass Points: 0 (for no) or 1 (for yes)
gravity_Yes = Lag_Struct_Params(14);          % Gravity: 0 (for no) or 1 (for yes)
%NOTE: Lag_Struct_Params(15),(16):            <- components of gravity vector (if gravity, initialize them below)
porous_Yes = Lag_Struct_Params(17);           % Porous Media: 0 (for no) or 1 (for yes)
electro_phys_Yes = Lag_Struct_Params(18);     % Electrophysiology (FitzHugh-Nagumo): 0 (for no) or 1 (for yes)
d_Springs_Yes = Lag_Struct_Params(19);        % Damped Springs: 0 (for no) or 1 (for yes)
update_D_Springs_Flag = Lag_Struct_Params(20);% Update_Damped_Springs: % 0 (for no) or 1 (for yes)
boussinesq_Yes = Lag_Struct_Params(21);       % Boussinesq Approx.: 0 (for no) or 1 (for yes)
exp_Coeff = Lag_Struct_Params(22);            % Expansion Coefficient (e.g., thermal, etc) for Boussinesq approx.
general_force_Yes = Lag_Struct_Params(23);    % General User-Defined Force Term: 0 (for no) or 1 (for yes)
poroelastic_Yes = Lag_Struct_Params(24);      % Poro-elastic Boundary: 0 (for no) or 1 (for yes)
coagulation_Yes = Lag_Struct_Params(25);      % Coagulation Model: 0 (for no) or 1 (for yes)
brinkman_Yes= Lag_Struct_Params(26);          % Brinkman fluid: 0 (for no) or 1 (for yes)

concentration_Yes = Con_Params{1};    % Background Concentration Gradient: 0 (for no) or 1 (for yes)
kDiff = Con_Params{2};                % Diffusion Coefficient
adv_scheme = Con_Params{3};           % Which advection scheme to use: 0 (for 1st O upwind) or 1 (for 3rd O WENO)
source_Yes = Con_Params{4};           % Do we have source model to use: 0 (for no) or 1 (for yes)
C_num = Con_Params{5};                % Number of Concentrations 

%--------------------------------------------
% Added in case no concentration being used;
%       sets concentration_Yes to 'none'
%--------------------------------------------
if isempty(concentration_Yes)
    concentration_Yes = 0;
end

% CLEAR INPUT DATA %
clear Fluid_Params Grid_Params Time_Params Lag_Name_Params Con_Params;


%Lagrangian Structure Data
ds = min( Lx/(2*Nx), Ly/(2*Ny) );   % Lagrangian Spacing
grid_Info(9) = ds;                  % Store Lagrangian resolution, ds


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
[idX,idY] = meshgrid(0:Nx-1,0:Ny-1);     % INITIALIZE FOR FLUID SOLVER FFT FUNCTION



% % % % % HOPEFULLY WHERE I CAN READ IN INFO!!! % % % % %




%--------------------------------------------------
% READ IN LAGRANGIAN POINTS %
%--------------------------------------------------
[Nb,xLag,yLag] = read_Vertex_Points(struct_name);
grid_Info(8) = Nb;          % # Total Number of Lagrangian Pts.
xLag_P = xLag;              % Initialize previous Lagrangian x-Values (for use in muscle-model)
yLag_P = yLag;              % Initialize previous Lagrangian y-Values (for use in muscle-model)


fprintf('\n--> FIBER MODEL INCLUDES: \n');


%--------------------------------------------------
% READ IN GEOMETRY CONNECTIONS (*if file exists*)
%--------------------------------------------------
test_ver = ver('MATLAB');
year_ver = test_ver.Release;
year = str2num(year_ver(3:6));
lett = year_ver(7);
if ((year<=2016 && strcmp(lett,'b'))||(( year<=2017 ) && ( strcmp(lett,'a')))) && ( exist([struct_name '.geo_connect'])== 2 ) 
    fprintf('  -Specified Geometric Connections\n');
    geo_Connect_MAT = read_Geometry_Connections(struct_name);
    flag_Geo_Connect = 1;
elseif isfile([struct_name '.geo_connect'])
    fprintf('  -Specified Geometric Connections\n');
    geo_Connect_MAT = read_Geometry_Connections(struct_name);
    flag_Geo_Connect = 1;
else
    geo_Connect_MAT = 0;
    flag_Geo_Connect = 0;
end
clear test_ver year_ver year lett;




%--------------------------------------------------
% READ IN SPRINGS (IF THERE ARE SPRINGS) %
%--------------------------------------------------
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





%-------------------------------------------------------------
% READ IN BEAMS (IF THERE ARE BEAMS aka TORSIONAL SPRINGS) %
%-------------------------------------------------------------
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





%----------------------------------------------------------------
% READ IN NON-INVARIANT BEAMS (IF THERE ARE NONINVARIANT BEAMS) %
%----------------------------------------------------------------
if ( nonInv_beams_Yes == 1)
    fprintf('  -Beams ("non-Invariant") and ');
    if update_nonInv_Beams_Flag == 0
        fprintf('NOT dynamically updating non-invariant beam properties\n');
    else
        fprintf('dynamically updating non-invariant beam properties\n');
    end
    
    nonInv_beams_info = read_nonInv_Beam_Points(struct_name);
    %beams:      col 1: 1ST PT.
    %            col 2: MIDDLE PT. (where force is exerted)
    %            col 3: 3RD PT.
    %            col 4: beam stiffness
    %            col 5: x-curavture
    %            col 6: y-curavture
else
    nonInv_beams_info = 0;
end








%---------------------------------------------------
% READ IN TARGET POINTS (IF THERE ARE TARGET PTS) %
%---------------------------------------------------
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



%--------------------------------------------------
% READ IN MASS POINTS (IF THERE ARE MASS PTS) %
%--------------------------------------------------
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

%-------------------------------------------------------
% CONSTRUCT GRAVITY INFORMATION (IF THERE IS GRAVITY) %
%-------------------------------------------------------
if gravity_Yes == 1
    %gravity_Vec(1) = Lag_Struct_Params(12);     % x-Component of Gravity Vector
    %gravity_Vec(2) = Lag_Struct_Params(13);     % y-Component of Gravity Vector
    xG = Lag_Struct_Params(15);
    yG = Lag_Struct_Params(16);
    normG = sqrt( xG^2 + yG^2 );
    gravity_Info = [gravity_Yes xG/normG yG/normG];
    %   col 1: flag if considering gravity
    %   col 2: x-component of gravity vector (normalized)
    %   col 3: y-component of gravity vector (normalized)
    
    clear xG yG normG;
    
else
    gravity_Info = 0;
end

%----------------------------------------------------
% READ IN POROUS MEDIA INFO (IF THERE IS POROSITY) %
%----------------------------------------------------
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




%----------------------------------------------------------------
% READ IN PORO-ELASTIC MEDIA INFO (IF THERE IS PORO-ELASTICITY) %
%----------------------------------------------------------------
if ( poroelastic_Yes == 1)
    fprintf('  -Poroelastic media\n');
    poroelastic_info = read_PoroElastic_Points(struct_name);
    F_Poro = zeros( length( poroelastic_info(:,1) ), 2);   % Initialization
    %poroelastic_info: col 1: Lag Pt. ID w/ Associated Porous Pt.
    %                  col 2: Porosity coefficient
else
    poroelastic_info = 0;
    F_Poro = 0;
end





%--------------------------------------------------
% READ IN DAMPED SPRINGS (IF THERE ARE SPRINGS) %
%--------------------------------------------------
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







%--------------------------------------------------
% READ IN MUSCLES (IF THERE ARE MUSCLES) %
%--------------------------------------------------
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











%--------------------------------------------------
% READ IN MUSCLES (IF THERE ARE MUSCLES) %
%--------------------------------------------------
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






%----------------------------------------------------------------------
% READ IN GENERAL FORCE PARAMETERS (IF THERE IS A USER-DEFINED FORCE) %
%----------------------------------------------------------------------
if ( general_force_Yes == 1 )
    fprintf('  -GENERAL FORCE MODEL (user-defined force term)\n');
    gen_force_info = read_General_Forcing_Function(struct_name);
        %
        %
        % ALL PARAMETERS / FORCE FUNCTION SET BY USER!
        %
        %
else
    gen_force_info = 0;  %just to pass placeholder into "please_Find_Lagrangian_Forces_On_Eulerian_grid function"
end





%-----------------------------------------------------------
% READ IN COAGULATION PARAMETERS (IF THERE IS COAGULATION) %
%-----------------------------------------------------------
if ( coagulation_Yes == 1 )
    fprintf('  -COAGULATION MODEL included\n');
    coagulation_info = read_Coagulation_Input_File(struct_name);
    aggregate_list = 0;
        %
        %
        % ALL PARAMETERS / FORCE FUNCTION SET BY USER!
        %
        %
else
    coagulation_info = 0;  %just to pass placeholder into "please_Find_Lagrangian_Forces_On_Eulerian_grid function"
    aggregate_list = 0;
end





%--------------------------------------------------
% SOLVE ELECTROPHYSIOLOGY MODEL FOR PUMPING %
%--------------------------------------------------
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

% PRINT IF PRESCRIBED BACKGROUND FLOW %
if ( arb_ext_force_Yes == 1)
fprintf('  -Prescribed background flow w/ penalty force\n');    
end

% PRINT IF TRACERS %
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

%------------------------------------------------------------------
% READ IN CONCENTRATIONS (IF THERE IS A(RE) BACKGROUND CONCENTRATION(S)) %
%------------------------------------------------------------------
if ( concentration_Yes == 1 )
    fprintf('  -Background concentration included\n');
    if C_num<2
    [C] = read_In_Concentration_Info(struct_name,Nx,0);
    else
        for ic=1:C_num
            [C1] = read_In_Concentration_Info(struct_name,Nx,ic);
            C(:,:,ic)=C1;
        end
    end
        %C:           Initial background concentration
        %kDiffusion:  Diffusion constant for Advection-Diffusion
else
    C = 0; % placeholder for plotting 
end


%--------------------------------------------------------------
% READ IN BRINKMAN PERMEABILITY (IF THERE IS A BRINKMAN TERM) %
%--------------------------------------------------------------
if ( brinkman_Yes == 1 )
    fprintf('  -Brinkman Term included\n');
    [Brink,kDiffusion] = read_In_Brinkman_Info(struct_name,Nx);
        %Brink:   Initial brinkman permeability matrix
        %kBrink:  Permeability for brinkman model
else
    Brink = 0; % placeholder for plotting 
end



%---------------------------------------------------------
% CONSTRUCT BOUSSINESQ INFORMATION (IF USING BOUSSINESQ) %
%---------------------------------------------------------
if boussinesq_Yes == 1
    fprintf('  -Boussinesq Approximation included\n');
    if exp_Coeff == 0
        fprintf('    -> exp_Coeff set to 1.0 by default, was assigned 0 in input2d <-\n');
        exp_Coeff = 1.0;
    end
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
    %Cxx = DD(C,dx,'x');
    %Cyy = DD(C,dy,'y');
    %laplacian_C = Cxx+Cyy;
    %C_0 = zeros(size(C)); % Define background concentration
end








% Initialize external fluid forcing flags
if arb_ext_force_Yes == 1 
    firstExtForce = 1;                           % initialize external forcing
    indsExtForce = 0;                            % initialize for external forcing computation
end

% ACTUAL TIME-STEPPING IBM SCHEME! 
%(flags for storing structure connects for printing and printing to .vtk)
cter = 0; ctsave = 0; firstPrint = 1; loc = 1; diffy = 1;


% CREATE VIZ_IB2D/HIER_IB2d_DATA FOLDER for .VTK FILES
mkdir('viz_IB2d');
if Output_Params(18) == 1
    mkdir('hier_IB2d_data'); 
end

% PRINT BRINKMAN PERMEABILITY TO MATRIX (*if BRINKMAN TERM*) %
if brinkman_Yes == 1
    cd('viz_IB2d');
    strNUM = give_String_Number_For_VTK(ctsave);
    brinkName = ['Brink.' strNUM '.vtk'];
    savevtk_scalar(Brink', brinkName, 'Brink',dx,dy);
    clear brinkName strNUM;
    cd ..;
end


%
% INITIALIZES SIMULATION -> if RESTARTing simulation, reads in previous
%                           data. Needs script to do so in Example/ folder
%
if restart_Flag == 0

    % Initialize the initial velocities to zero.
    U = zeros(Ny,Nx);                                % x-Eulerian grid velocity
    V = U;                                           % y-Eulerian grid velocity
    mVelocity = zeros( length(mass_info(:,1)), 2 );  % mass-Pt velocity 

    %Initialize Vorticity, uMagnitude, and Pressure for initial colormap
    %Print initializations to .vtk
    vort=zeros(Nx,Ny); uMag=vort; p = vort;  lagPts = [xLag yLag zeros(Nb,1)]; 
    [connectsMat,spacing] = give_Me_Lag_Pt_Connects(ds,xLag,yLag,Nx,springs_Yes,springs_info);
    [dconnectsMat,dspacing] = give_Me_Lag_Pt_Connects(ds,xLag,yLag,Nx,d_Springs_Yes,d_springs_info);
    Fxh = vort; Fyh = vort; F_Lag = zeros( Nb, 2); 
    print_vtk_files(Output_Params,0,ctsave,vort,uMag,p,U',V',Lx,Ly,Nx,Ny,lagPts,springs_Yes,connectsMat,tracers,concentration_Yes,C,Fxh,Fyh,F_Lag,coagulation_Yes,aggregate_list,d_Springs_Yes,dconnectsMat,poroelastic_Yes,F_Poro,C_num);
    fprintf('\n |****** Begin IMMERSED BOUNDARY SIMULATION! ******| \n\n');
    fprintf('Current Time(s): %6.6f\n\n',current_time);
    ctsave = ctsave+1;

else
    
    %
    % PASSES BACK PREVIOUS: 1. uVec components (U,V)
    %                       2. last xLag and yLag 
    %                       3. 2nd to last xLag and yLag
    %                       4. sets counters, current_time, saving counters
    %
    [current_time,cter,ctsave,U,V,xLag,yLag,xLag_P,yLag_P,path_to_Viz] = help_Me_Restart(dt);
    clear restart_Flag;
    
    %
    % RESTART CONCENTRATION 
    %       
    %       note: need ctsave - 1, b/c it was incremented in help_Me_Restart()
    %
    if concentration_Yes == 1
        C = pass_Back_Concentration_Data_For_Restart(ctsave-1,path_to_Viz);
    else 
        C = 0;
    end
    
    
    %
    % RESTART MASS * STILL NEEDS IMPLEMENTATION <-> default set *
    %
    mVelocity = zeros( length(mass_info(:,1)), 2 );  % mass-Pt velocity  
    
    %
    % RESTART POROELASTIC * STILL NEEDS IMPLEMENTATION <-> default set *
    %
    if poroelastic_Yes == 1
        F_Poro = pass_Back_Data_To_Restart(path_to_Viz,'poroelastic',ctsave-1);
    else 
        F_Poro = 0;
    end
    
    %[~, ~, ~, ~, F_Poro, ~] = please_Find_Lagrangian_Forces_On_Eulerian_grid(dt, current_time, xLag_h, yLag_h, xLag_P, yLag_P, x, y, grid_Info, Lag_Struct_Params, springs_info, target_info, beams_info, nonInv_beams_info ,muscles_info, muscles3_info, mass_info, electro_potential, d_springs_info, gen_force_info, poroelastic_info, coagulation_info, aggregate_list);

    
    %
    % INITIALIZATION
    %
    vort=zeros(Nx,Ny); uMag=vort; p = vort;  lagPts = [xLag yLag zeros(Nb,1)]; 
    [connectsMat,spacing] = give_Me_Lag_Pt_Connects(ds,xLag,yLag,Nx,springs_Yes,springs_info);
    [dconnectsMat,dspacing] = give_Me_Lag_Pt_Connects(ds,xLag,yLag,Nx,d_Springs_Yes,d_springs_info);
    Fxh = vort; Fyh = vort; F_Lag = zeros( Nb, 2); 
    
    fprintf('\n |****** RESTARTED IMMERSED BOUNDARY SIMULATION! ******| \n\n');
    fprintf('Start Time(s): %6.6f\n\n',current_time-dt); % - dt just to keep with previous printing convention
    
end


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
   
    %--------------------------------------------------------------------
    % DEFINE VARIABLES TO HOLD PREVIOUS DATA FOR CONCENTRATION DYNAMICS
    %--------------------------------------------------------------------
    if concentration_Yes == 1 && source_Yes(1)>0
        xLag_prev=xLag;
        yLag_prev=yLag;
        U_prev=U;
        V_prev=V;
    elseif concentration_Yes == 1
        U_prev=U;
        V_prev=V;
    end
    
    %--------------------------------------------------------------------------------
    %
    %**************** Step 1: Update Position of Boundary of membrane at half time-step *******************
    %                           (Variables end with h if it is a half-step)
    %
    %--------------------------------------------------------------------------------
    [xLag_h, yLag_h] = please_Move_Lagrangian_Point_Positions(mu, U, V, xLag, yLag, xLag, yLag, x, y, dt/2, grid_Info,0,poroelastic_Yes,poroelastic_info,F_Poro);
    
    if mass_Yes == 1
       [mass_info, massLagsOld] = please_Move_Massive_Boundary(dt/2,mass_info,mVelocity); 
    end
    
    if ( ( electro_phys_Yes == 1) && (muscles_Yes == 0) ) %0.75e0 %1e1 for Re=3.4    %2.25e1 for Re=0.5 D=10
        springs_info(ePhys_Start:ePhys_End,3) = ( 1.6e1*electro_potential(ePhys_Ct,:)') .^4;%( 1e4*electro_potential(ePhys_Ct,:)'.*springs_info(ePhys_Start:ePhys_End,3) ).^4;
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
    
    if ( ( update_nonInv_Beams_Flag == 1 ) && ( nonInv_beams_Yes == 1) )
       nonInv_beams_info = update_nonInv_Beams(dt,current_time,nonInv_beams_info); 
    end
    
    if ( ( update_D_Springs_Flag == 1 ) && ( d_Springs_Yes == 1) )
       d_springs_info = update_Damped_Springs(dt,current_time,d_springs_info); 
    end
    
    
    %--------------------------------------------------------------------------------
    %
    %**************** STEP 2: Calculate Force coming from membrane at half time-step ****************
    %
    %--------------------------------------------------------------------------------
   [Fxh, Fyh, F_Mass_Bnd, F_Lag, F_Poro, aggregate_list] =    please_Find_Lagrangian_Forces_On_Eulerian_grid(dt, current_time, xLag_h, yLag_h, xLag_P, yLag_P, x, y, grid_Info, Lag_Struct_Params, springs_info, target_info, beams_info, nonInv_beams_info ,muscles_info, muscles3_info, mass_info, electro_potential, d_springs_info, gen_force_info, poroelastic_info, coagulation_info, aggregate_list, flag_Geo_Connect, geo_Connect_MAT);
    
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
    
    
    %--------------------------------------------------------------------------------
    %
    %***************** STEP 3: Solve for Fluid motion ******************************************
    %
    %--------------------------------------------------------------------------------
    
    % Add in effect from BOUSSINESQ
    if boussinesq_Yes == 1
        Fxh = Fxh + rho*fBouss_X.*(C)*(tanh(20*current_time));
        Fyh = Fyh + rho*fBouss_Y.*(C)*(tanh(20*current_time));
    end
    
    
    % Add in effect from BRINKMAN 
    if brinkman_Yes == 1
        Fxh = Fxh - mu*Brink.*U; % ADD BRINKMAN TERM!!!
        Fyh = Fyh - mu*Brink.*V; % ADD BRINKMAN TERM!!!
    end

    % Update fluid motion
    [Uh, Vh, U, V, p] =   please_Update_Fluid_Velocity(U, V, Fxh, Fyh, rho, mu, grid_Info, dt, idX, idY);
    
    
    %--------------------------------------------------------------------------------
    %
    %***************** STEP 4: Update Position of Boundary of membrane again for a half time-step *******
    %
    %--------------------------------------------------------------------------------
    xLag_P = xLag_h;   % Stores old Lagrangian x-Values (for muscle model)
    yLag_P = yLag_h;   % Stores old Lagrangian y-Values (for muscle model)
    %Uh, Vh instead of U,V?
    [xLag, yLag] =     please_Move_Lagrangian_Point_Positions(mu, Uh, Vh, xLag, yLag, xLag_h, yLag_h, x, y, dt, grid_Info,porous_Yes,poroelastic_Yes,poroelastic_info,F_Poro);
    

    %IF POROUS NODES, NOTE: SET UP FOR BOTH CLOSED + OPEN SYSTEMS NOW!!!
    if porous_Yes == 1
       [Por_Mat,nX,nY] = please_Compute_Porous_Slip_Velocity(ds,xLag,yLag,porous_info,F_Lag);
       xLag( porous_info(:,1) ) = xLag( porous_info(:,1) ) - dt*( Por_Mat(:,1)+Por_Mat(:,2) ).*nX;
       yLag( porous_info(:,1) ) = yLag( porous_info(:,1) ) - dt*( Por_Mat(:,1)+Por_Mat(:,2) ).*nY;
       porous_info(:,2) = xLag( porous_info(:,1) );    %Stores x-Lags as x-Porous Pt. Identities
       porous_info(:,3) = yLag( porous_info(:,1) );    %Stores y-Lags as y-Porous Pt. Identities
       xLag = mod(xLag, Lx);                           % If structure goes outside domain
       yLag = mod(yLag, Ly);                           % If structure goes outside domain
    end    
    
    
    % If there are tracers, update tracer positions %
    if tracers_Yes == 1
        %Uh, Vh instead of U,V?
        [xT, yT] = please_Move_Lagrangian_Point_Positions(mu, Uh, Vh, xT, yT, xT, yT, x, y, dt, grid_Info,0,0,0,0); %0 for always no porous tracers / poroelastic elements
        tracers(:,2) = xT;
        tracers(:,3) = yT;
    end
    
    % If there is a background concentration, update concentration-gradient %
    if concentration_Yes == 1 && source_Yes(1)==0
      
        % ** Advection-diffusion WITHOUT a source term **

        % Not fully integrated nor tested.
        %[C,~] = please_Update_Adv_Diff_Concentration_Flux_Limiter_FV(C,dt,dx,dy,U,V,kDiffusion); 
        
        % Update concentration (solve advection-diffusion eq) w/ either UPWIND
        %   or WENO advection scheme
        [N1,N2,num_con]=size(C);
        
        % iterate over concentrations
        for ic=1:num_con
            [C1,~] = please_Update_Adv_Diff_Concentration(C(:,:,ic),dt,dx,dy,U_prev,V_prev,kDiff(ic),adv_scheme,Lx,Ly);
            C(:,:,ic)=C1;
        end
        
    elseif concentration_Yes == 1 && source_Yes(1)>0
        
        % initialize info
        [N1,N2,num_con]=size(C);
        
        % interate over concentrations
        for ic=1:num_con
        % Compute force term for advection-diffusion with a source term
        
            fs = please_Find_Source_For_Concentration(dt, current_time, xLag_prev, yLag_prev, x, y, grid_Info, C, flag_Geo_Connect, geo_Connect_MAT,ic);

            % Update concentration (solver advection-diffusion eq) w/ either UPWIND
            %   or WENO advection scheme
            [C1,~] = please_Update_Adv_Diff_Concentration_Source(C(:,:,ic),dt,dx,dy,U_prev,V_prev,kDiff(ic),fs,Lx,Ly,adv_scheme);
            C(:,:,ic)=C1;
            
        end
        
   end
        
    % Plot Lagrangian/Eulerian Dynamics! %
    if ( ( mod(cter,pDump) == 0  ) && ( cter >= pDump ) )
        
        %Compute vorticity, uMagnitude
        vort = give_Me_Vorticity(U,V,dx,dy);
        uMag = give_Me_Magnitude_Velocity(U,V);
        
        %Plot in Matlab
        if pMatlab == 1
            [loc, diffy] = please_Plot_Results(ds,X,Y,U,V,vort,uMag,p,C,xLag,yLag,lagPlot,velPlot,vortPlot,pressPlot,uMagPlot,conPlot,firstPrint,loc,diffy,spacing);
        end
        
        %Print .vtk files!
        lagPts = [xLag yLag zeros(length(xLag),1)];
        print_vtk_files(Output_Params,current_time,ctsave,vort,uMag',p',U',V',Lx,Ly,Nx,Ny,lagPts,springs_Yes,connectsMat,tracers,concentration_Yes,C,Fxh',Fyh',F_Lag,coagulation_Yes,aggregate_list,d_Springs_Yes,dconnectsMat,poroelastic_Yes,F_Poro,C_num);
        
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
        
    end

    
    % Update current_time value & counter
    current_time = current_time+dt;
    cter = cter + 1;
    
end %ENDS TIME-STEPPING LOOP










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives appropriate string number for filename in printing the
% .vtk files. for viz_IB2d and hier_IB2d_data folders
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_vtk_files(Output_Params,current_time,ctsave,vort,uMag,p,U,V,Lx,Ly,Nx,Ny,lagPts,springs_Yes,connectsMat,tracers,concentration_Yes,C,fXGrid,fYGrid,F_Lag,coag_Yes,aggregate_list,d_springs_Yes,dconnectsMat,poroelastic_Yes,F_Poro,C_n)

    %
    %  Output_Params(1):  print_dump
    %               (2):  plot_Matlab
    %               (3):  plot_LagPts
    %               (4):  plot_Velocity
    %               (5):  plot_Vorticity
    %               (6):  plot_MagVelocity
    %               (7):  plot_Pressure
    %               (8):  save_Vorticity 
    %               (9):  save_Pressure 
    %               (10): save_uVec 
    %               (11): save_uMag 
    %               (12): save_uX 
    %               (13): save_uY 
    %               (14): save_fMag 
    %               (15): save_fX 
    %               (16): save_fY 
    %               (17): save_hier 

    
%Give EULERIAN spacing for grid
dx = Lx/Nx; 
dy = Ly/Ny;

%
% Print all Eulerian Data and Lagrangian Position Data to viz_IB2d folder
%
cd('viz_IB2d'); %Go into viz_IB2d directory

    %Find string number for storing files
    strNUM = give_String_Number_For_VTK(ctsave);
    lagPtsName = ['lagsPts.' strNUM '.vtk'];

    %Print Lagrangian Pts to .vtk format
    savevtk_points(lagPts, lagPtsName, 'lagPts', current_time);

    % Print Spring Connections (* if springs *)
    if springs_Yes == 1
        %Print Lagrangian Pts w/ CONNECTIONS to .vtk format
        lagPtsConName=['lagPtsConnect.' strNUM '.vtk'];
        %connectsMat2 = test_Connects_Distances(lagPts(:,1),lagPts(:,2),connectsMat);
        connectsMat2 = connectsMat;
        savevtk_points_connects(lagPts, lagPtsConName, 'lagPtsConnected',connectsMat2, current_time);
    end
    
    % Print DAMP Spring Connections (* if damped springs *)
    if d_springs_Yes == 1 
        %Print Lagrangian Pts w/ CONNECTIONS to .vtk format
        dlagPtsConName=['damp_Connect.' strNUM '.vtk'];
        %dconnectsMat2 = test_Connects_Distances(lagPts(:,1),lagPts(:,2),dconnectsMat);
        dconnectsMat2 = dconnectsMat;
        savevtk_points_connects(lagPts, dlagPtsConName, 'dampConnected',dconnectsMat2, current_time);
    end
    
    %if (  ( coag_Yes == 1 ) && ( aggregate_list(1) > 0 ) )
    if ( coag_Yes == 1 )
        %Print COAGULATION CONNECTIONS to .vtk format
        CoagConName=['Coag_Connect.' strNUM '.vtk'];
        if aggregate_list(1) > 0
            agg_list_aux = test_Connects_Distances(lagPts(:,1),lagPts(:,2),aggregate_list(:,3:4));
            savevtk_points_connects(lagPts, CoagConName, 'CoagConnected',[connectsMat2; agg_list_aux], current_time);
        else
            savevtk_points_connects(lagPts, CoagConName, 'CoagConnected',connectsMat2, current_time);
        end
    end
    
    clear connectsMat2;
    
    %Print Tracer Pts (*if tracers*)
    if tracers(1,1) == 1
        tracersPtsName = ['tracer.' strNUM '.vtk'];
        savevtk_points(tracers(:,2:4),tracersPtsName, 'tracers', current_time); 
    end


    %Print SCALAR DATA (i.e., colormap data) to .vtk file
    if Output_Params(9) == 1
        vortfName = ['Omega.' strNUM '.vtk'];
        savevtk_scalar(vort, vortfName, 'Omega',dx,dy, current_time);
    end
    if Output_Params(10) == 1
        pfName = ['P.' strNUM '.vtk'];
        savevtk_scalar(p, pfName, 'P',dx,dy, current_time);
    end
    if Output_Params(12) == 1
        uMagfName = ['uMag.' strNUM '.vtk'];
        savevtk_scalar(uMag, uMagfName, 'uMag',dx,dy, current_time);
    end
    if Output_Params(13) == 1
        uXName = ['uX.' strNUM '.vtk'];
        savevtk_scalar(U, uXName, 'uX',dx,dy, current_time);
    end
    if Output_Params(14) == 1
        uYName = ['uY.' strNUM '.vtk'];
        savevtk_scalar(V, uYName, 'uY',dx,dy, current_time);
    end
    if Output_Params(16) == 1
        fXName = ['fX.' strNUM '.vtk'];
        savevtk_scalar(fXGrid, fXName, 'fX',dx,dy, current_time);
    end
    if Output_Params(17) == 1
        fYName = ['fY.' strNUM '.vtk'];
        savevtk_scalar(fYGrid, fYName, 'fY',dx,dy, current_time);
    end
    if Output_Params(15) == 1
        fMagName = ['fMag.' strNUM '.vtk'];
        savevtk_scalar(sqrt( fXGrid.^2 + fYGrid.^2 ), fMagName, 'fMag',dx,dy, current_time);
    end
    
    %
    % FOR RESTARTING SIMULATION 
    %
    % If POROELASTIC Structure (saves last F_Poro for restart functionality) %
    if poroelastic_Yes == 1
    
        if ctsave == 0
            mkdir('Data_For_Restart');
        end
        
        cd('Data_For_Restart');
        print_Matrix_to_TXT_for_Restart(['poroelastic_Restart_' num2str(ctsave) '.txt'],F_Poro);
        cd ../; % Return to viz_IB2d folder
        
    end
    
    
    % If Background Concentration Gradient %
    if concentration_Yes == 1
        %interate over concentrations
        
        if C_n==0               
            confName = [sprintf('concentration.') strNUM '.vtk'];
            inds=C<1e-15; % Gives logical matrix
            C(inds) = 0;  % Zeros out values of matrix, C
            savevtk_scalar(C', confName, 'Concentration',dx,dy, current_time);
        else
            for ic=1:C_n
                confName = [sprintf('concentration%g.',ic) strNUM '.vtk'];
                inds=C<1e-15; % Gives logical matrix
                C(inds) = 0;  % Zeros out values of matrix, C
                savevtk_scalar(C(:,:,ic)', confName, 'Concentration',dx,dy, current_time);
            end
        end
    end


    %Print VECTOR DATA (i.e., velocity data) to .vtk file
    if Output_Params(11) == 1
        velocityName = ['u.' strNUM '.vtk'];
        savevtk_vector(U, V, velocityName, 'u',dx,dy, current_time);
    end

%Get out of viz_IB2d folder
cd ..

if Output_Params(18) == 1
    if length( lagPts ) <= 5
        %
        % Print Lagrangian Force Data to hier_IB2d_data folder (if <= 5 lag pts)
        %
        cd('hier_IB2d_data'); %change directory to hier-data folder

            % Save x-y force data!
            fLag_XName = ['fX_Lag.' strNUM '.vtk'];
            fLag_YName = ['fY_Lag.' strNUM '.vtk'];
            savevtk_points_with_scalar_data( lagPts, F_Lag(:,1), fLag_XName, 'fX_Lag', current_time);
            savevtk_points_with_scalar_data( lagPts, F_Lag(:,2), fLag_YName, 'fY_Lag', current_time);

            % Save force Magnitude (no normal/tangential -> not enough points)
            fMagName = ['fMag.' strNUM '.vtk'];
            fLagMag = sqrt( F_Lag(:,1).^2 + F_Lag(:,2).^2 ); % Compute magnitude of forces on boundary
            savevtk_points_with_scalar_data( lagPts, fLagMag, fMagName, 'fMag', current_time);
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
            savevtk_points_with_scalar_data( lagPts, F_Lag(:,1), fLag_XName, 'fX_Lag', current_time);
            savevtk_points_with_scalar_data( lagPts, F_Lag(:,2), fLag_YName, 'fY_Lag', current_time);

            % Save force Magnitude, Mag. Normal, Mag. Tangential
            fMagName = ['fMag.' strNUM '.vtk'];
            fNormalName = ['fNorm.' strNUM '.vtk'];
            fTangentName = ['fTan.' strNUM '.vtk'];

            fLagMag = sqrt( F_Lag(:,1).^2 + F_Lag(:,2).^2 ); % Compute magnitude of forces on boundary

            savevtk_points_with_scalar_data( lagPts, fLagMag, fMagName, 'fMag', current_time);
            savevtk_points_with_scalar_data( lagPts, F_Normal_Mag, fNormalName, 'fNorm', current_time);
            savevtk_points_with_scalar_data( lagPts, F_Tan_Mag, fTangentName, 'fTan', current_time);

        cd .. % Get out of hier_IB2d_data folder
    end
end % ENDS IF-STATEMENT FOR IF SAVE_HIER

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: returns appropriate (x,y) values to analyze inside the
% ventricle
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Matrix_to_TXT_for_Restart(strTitle,matrix)

fid = fopen(strTitle, 'wt'); % Open for writing
fprintf(fid, 'poroelastic\n');
for i=1:size(matrix,1)
   fprintf(fid, '%d ', matrix(i,:));
   fprintf(fid, '\n');
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: tests ConnectsMat for too far of distances to print connection!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function connectsMat2 = test_Connects_Distances(xLag,yLag,connectsMat)

Lx = 1;
r = 0.2*Lx/18.1;
connectsMat = connectsMat + 1;
j=1;

for i=1:length( connectsMat(:,1) )
    ind1 = connectsMat(i,1);
    ind2 = connectsMat(i,2);
    dist = sqrt( ( xLag(ind1)-xLag(ind2) )^2 + ( yLag(ind1)- yLag(ind2) )^2   );
    if dist < 5*r
        connectsMat2(j,1) = ind1-1;  % minus for Cpp notation
        connectsMat2(j,2) = ind2-1;  % minus for Cpp notation
        j=j+1;
    end
end



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

function savevtk_points_connects( X, filename, vectorName,connectsMat, time)

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
fprintf(file, 'FIELD FieldData 1\n');
fprintf(file, 'TIME 1 1 double\n');
fprintf(file, '%.8f\n',time);
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

function savevtk_points( X, filename, vectorName, time)

%X is matrix of size Nx3

N = length( X(:,1) );


%TRY PRINTING THEM AS UNSTRUCTURED_GRID
file = fopen (filename, 'w');
fprintf(file, '# vtk DataFile Version 2.0\n');
fprintf(file, [vectorName '\n']);
fprintf(file, 'ASCII\n');
fprintf(file, 'DATASET UNSTRUCTURED_GRID\n\n');
%
fprintf(file, 'FIELD FieldData 1\n');
fprintf(file, 'TIME 1 1 double\n');
fprintf(file, '%.8f\n',time);
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

function savevtk_vector(X, Y, filename, vectorName,dx,dy, time)
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
    %
    fprintf(fid, 'FIELD FieldData 1\n');
    fprintf(fid, 'TIME 1 1 double\n');
    fprintf(fid, '%.8f\n',time);
    %
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

function savevtk_scalar(array, filename, colorMap,dx,dy, time)
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
    %
    fprintf(fid, 'FIELD FieldData 1\n');
    fprintf(fid, 'TIME 1 1 double\n');
    fprintf(fid, '%.8f\n',time);
    %
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

function savevtk_points_with_scalar_data( X, scalarArray, filename, vectorName, time)

%X is matrix of size Nx3

N = length( X(:,1) );


%TRY PRINTING THEM AS UNSTRUCTURED_GRID
file = fopen (filename, 'w');
fprintf(file, '# vtk DataFile Version 2.0\n');
fprintf(file, [vectorName '\n']);
fprintf(file, 'ASCII\n');
fprintf(file, 'DATASET UNSTRUCTURED_GRID\n\n');
%
fprintf(file, 'FIELD FieldData 1\n');
fprintf(file, 'TIME 1 1 double\n');
fprintf(file, '%.8f\n',time);
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
% FUNCTION: Reads in Geometry Connections Between Lag IDs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function geo_MAT = read_Geometry_Connections(struct_name)

filename = [struct_name '.geo_connect'];  %Name of file to read in
fileID = fopen(filename);

% Read in the file, use 'CollectOutput' to gather all similar data together
% and 'CommentStyle' to to end and be able to skip lines in file.
C = textscan(fileID,'%f %f','CollectOutput',1);

fclose(fileID);     %Close the data file.

geo_MAT = C{1};    % Stores all read in data





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
%           (also used to read in permeability matrix for Brinkman &
%           permeability Brinkman constant, k.)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C] = read_In_Concentration_Info(struct_name,N,nc)
if nc==0
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
    %kDiff = con_info(1,1);     %coefficient of diffusion
    C = con_info;%(2:end,1:end); %initial concentration
else
    % reads in multiple concentrations
    filename = [struct_name sprintf('.concentration_%g',nc)];  %Name of file to read in

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
    %kDiff = con_info(1,1);     %coefficient of diffusion
    C = con_info;%(2:end,1:end); %initial concentration

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in the Brinkman permeability and initial permeability,
% Brinkman
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Brinkman,kBrink] = read_In_Brinkman_Info(struct_name,N)

filename = [struct_name '.brinkman'];  %Name of file to read in

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
kBrink = con_info(1,1);     %coefficient of diffusion
Brinkman = con_info(2:end,1:end); %initial concentration






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
%         col 3: flag for stencil point







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in the # of POROELASTIC PTS, PORO-ELASTIC PT-NODEs, and 
%           their POROUSITY-COEFFICIENTS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function porosity = read_PoroElastic_Points(struct_name)

filename = [struct_name '.poroelastic'];  %Name of file to read in

fileID = fopen(filename);

    % Read in the file, use 'CollectOutput' to gather all similar data together
    % and 'CommentStyle' to to end and be able to skip lines in file.
    C = textscan(fileID,'%f %f','CollectOutput',1);

fclose(fileID);        %Close the data file.

porous_info = C{1};    %Stores all read in data in vertices (N+1,2) array

%Store all elements on .spring file into a matrix starting w/ 2nd row of read in data.
porosity = porous_info(2:end,1:2);

%poroelastic:  col 1: Lag Pt. ID w/ Associated Porous Pt.
%              col 2: Porosity coefficient










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
% FUNCTION: Reads in the # of NON-INVARIANT beams and all 1st Pt, MIDDLE Pt, 
%           and 3rd Pt beam STIFFNESSES, and CURVATURE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function beams = read_nonInv_Beam_Points(struct_name)

filename = [struct_name '.nonInv_beam'];  %Name of file to read in

fileID = fopen(filename);

    % Read in the file, use 'CollectOutput' to gather all similar data together
    % and 'CommentStyle' to to end and be able to skip lines in file.
    C = textscan(fileID,'%f %f %f %f %f %f','CollectOutput',1);

fclose(fileID);      %Close the data file.

beam_info = C{1};    %Stores all read in data in vertices (N+1,2) array

%Store all elements on .beam file into a matrix starting w/ 2nd row of read in data.
beams = beam_info(2:end,1:6);

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: READS IN ALL THE DATA FOR THE USER-DEFINED FORCE FUNCTION!!!
%           NOTE: DOES NOT SPECIFY HOW MANY PARAMETERS THERE ARE.
%           NOTE: COMPLETELY USER-DEFINED
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function force_general = read_General_Forcing_Function(struct_name)

filename = [struct_name '.user_force'];  %Name of file to read in

% Imports all the data into a data structure
gen_force_info = importdata(filename,' ',1);

% Save data to new variable, force_general
force_general = gen_force_info.data;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: returns matrix of data for restart of fiber model (e.g.,
%           poroelastic or mass, etc.)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MAT = pass_Back_Data_To_Restart(path_to_Viz,strName,ctsave)

% Grab correct data
filename = [strName '_Restart_' num2str(ctsave) '.txt'];

% Add path to viz_IB2d folder desired
addpath([path_to_Viz '/Data_For_Restart/']);

% Imports all the data into a data structure
matrix_info = importdata(filename,' ',1);

% Save data to new variable, MAT
MAT = matrix_info.data;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: READS IN ALL THE DATA FOR THE COAGULATION MODEL
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function coag = read_Coagulation_Input_File(struct_name)

filename = [struct_name '.coagulation'];  %Name of file to read in


fileID = fopen(filename);

    % Read in the file, use 'CollectOutput' to gather all similar data together
    % and 'CommentStyle' to to end and be able to skip lines in file.
    C = textscan(fileID,'%f','CollectOutput',1);

fclose(fileID);      % Close the data file.

coag_Info = C{1};    % Stores all read in data

%Store all elements on .coagulation file into a matrix starting w/ 2nd row of read in data.
coag = coag_Info(1:end,1);

% coag:     row 1: staring index for coagulation cells
%           row 2: threshold radii
%           row 3: threshold force
%           row 4: row N+1: # of Lag. Pts in Cell, 1
%           row 5: row N+1: # of Lag. Pts in Cell, 2
%            .  
%            .  
%            . 
%           row N+3: # of Lag. Pts in Cell, N

    
