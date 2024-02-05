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
% 	2. Beams (*torsional springs)
% 	3. Target Points
% 	4. Muscle-Model (combined Force-Length-Velocity model, "Hill+(Length-Tension)")
%   .
%   .
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Computes the components of the force term in Navier-Stokes from
%           deformations of the boundary of the immersed boundary
%
%       NOTE: (1) Commented out implementations illustrating attempts to make code more efficient
%             (2) Lots of old implementation included (but commented out) for teaching purposes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Fx, Fy, F_Mass, F_Lag, F_Poro, aggregate_list] = please_Find_Lagrangian_Forces_On_Eulerian_grid(dt, current_time, xLag, yLag,xLag_P,yLag_P, x, y, grid_Info, model_Info, springs, targets, beams, nonInv_beams, muscles, muscles3, masses, electro_potential, d_Springs, general_force,poroelastic_info, coagulation, aggregate_list, flag_Geo_Connect, geo_Connect_MAT)


%-----------------------------------------------------------------------------------------
% The components of the force are given by
%       F(x,y) = int{ f(s) * delta(x - xLag(s)) * delta(y - yLag(s)) * ds }
%                   where s parameteriizes the Lagrangian structure.
%
% xLag:           x positions of Lagrangian structure
% yLag:           y positions of Lagrangian structure  
% x:              x positions on Eulerian grid
% y:              y positions on Eulerian grid
% grid_Info:      holds lots of geometric pieces about grid / simulations
% model_Info:     Stores if springs, if update_springs, if target_pts, if update_target_pts (as 0 (no) or 1 (yes) )
% springs:        Stores Leader Node, Follower Node, Spring Stiffness, Restling-Lengths, all in column vecs
% beams:          Stores 1st Node, 2nd (MIDDLE-MAIN) Node, 3rd Nodes, Beam Stiffnesses, and Beam curvatures
% targets:        Stores target point index, correponding xLag, yLag and target point stiffness
% masses:         Stores mass point index, correponding xLag, yLag, "spring" stiffness, and mass value parameter
%    .
%    .
% coagulation:      Stores coagulation data: first index of lag pt. of cell, threshold radii, fracture force, # of pts in each cell
% aggregate_list:   Stores list of bonds between Lag. Pt. Indices (INDICES OF CELLS AS A WHOLE)
% flag_Geo_Connect: If User provides geometrical details about neighboring pts
% geo_Connect_MAT:  Gives which LAG IDs are geometrical neighbors
% current_time:     Current time of simulation (in seconds)
%
%-----------------------------------------------------------------------------------------

% Grid Info %
Nx =    grid_Info(1); % # of Eulerian pts. in x-direction
Ny =    grid_Info(2); % # of Eulerian pts. in y-direction
Lx =    grid_Info(3); % Length of Eulerian grid in x-coordinate
Ly =    grid_Info(4); % Length of Eulerian grid in y-coordinate
dx =    grid_Info(5); % Spatial-size in x
dy =    grid_Info(6); % Spatial-size in y
supp =  grid_Info(7); % Delta-function support
Nb =    grid_Info(8); % # of Lagrangian pts. 
ds =    grid_Info(9); % Lagrangian spacing


% Model Potential Forces %
springs_Yes = model_Info(1);        % Springs: 0 (for no) or 1 (for yes) 
target_pts_Yes = model_Info(3);     % Target_Pts: 0 (for no) or 1 (for yes)
beams_Yes = model_Info(5);          % Beams (Torsional Springs): 0 (for no) or 1 (for yes)
nonInv_beams_Yes = model_Info(7);
muscle_LT_FV_Yes = model_Info(9);   % Length-Tension/Force-Velocity Muscle: 0 (for no) or 1 (for yes) (Length/Tension - Hill Model)
muscle_3_Hill_Yes = model_Info(10); % 3-Element Hill Model: 0 (for no) or 1 (for yes) (3 Element Hill + Length-Tension/Force-Velocity)
mass_Yes = model_Info(13);          % Mass Pts: 0 (for no) or 1 (for yes)
electro_phys_Yes = model_Info(18);  % Electrophysiology (FitzHugh-Nagumo): 0 (for no) or 1 (for yes)
d_Springs_Yes = model_Info(19);     % Damped Springs: 0 (for no) or 1 (for yes)
gen_force_Yes = model_Info(23);     % General User-Defined Force: 0 (for no) or 1 (for yes)
poroelastic_Yes = model_Info(24);   % Poroelastic Media: 0 (for no) or 1 (for yes)
coagulation_Yes = model_Info(25);   % Coagulation Model: 0 (for no) or 1 (for yes)


%----------------------------------------------------------------------------------------------
% Compute MUSCLE LENGTH-TENSION/FORCE-VELOCITY (if using combined length/tension-Hill model) %
%----------------------------------------------------------------------------------------------
%
if ( ( muscle_LT_FV_Yes == 1 ) && ( electro_phys_Yes == 0 ) )
    [fx_muscles, fy_muscles] = give_Muscle_Force_Densities(Nb,xLag,yLag,xLag_P,yLag_P,muscles,current_time,dt);
elseif ( ( muscle_LT_FV_Yes == 1 ) && ( electro_phys_Yes == 1 ) )
    [fx_muscles, fy_muscles] = give_ElectroPhys_Ca_Muscle_Force_Densities(Nb,xLag,yLag,xLag_P,yLag_P,muscles,current_time,dt,electro_potential);
else
    fx_muscles = 0;%zeros(length(xLag),1); % No x-Forces coming from muscles
    fy_muscles = fx_muscles;               % No y-Forces coming from muscles
end



%-----------------------------------------------------------------------------------------
% Compute 3-ELEMENT HILL MUSCLE MODEL FORCE DENSITIES (if using combined 3-HILL + LT/FV) %
%-----------------------------------------------------------------------------------------
%
if ( muscle_3_Hill_Yes == 1)
    [fx_muscles3, fy_muscles3] = give_3_Element_Muscle_Force_Densities(Nb,xLag,yLag,xLag_P,yLag_P,muscles3,current_time,dt);
else
    fx_muscles3 = 0;%zeros(length(xLag),1); % No x-Forces coming from muscles
    fy_muscles3 = fx_muscles3;              % No y-Forces coming from muscles
end



%----------------------------------------------------------------------
% Compute SPRING FORCE DENSITIES (if there are springs!)
%----------------------------------------------------------------------
if ( springs_Yes == 1 )
    
    % Compute "Connections Matrix" for what springs are attached to whom     %
    % (and "Connections Stiffness Matrix" and "Connections Restling Lengths" %
    %connects = give_Me_Spring_Connections_Matrix(Nb,Nsprings,sp_1,sp_2,K_Vec,L_Vec); 
    
    % Compute distances between Lag-Pts w/ Springs for Spring-Tension Calc.
    %dLag_x = give_Spring_Lagrangian_Distance(xLag, Lx, springs);
    %dLag_y = give_Spring_Lagrangian_Distance(yLag, Ly, springs);

    % Compute the Lagrangian SPRING tensions!
    %[Tx Ty] = give_Me_Spring_Lagrangian_Tension(Nb,dLag_x,dLag_y,springs);

    % Compute the Lagrangian SPRING force densities!
    [fx_springs, fy_springs] = give_Me_Spring_Lagrangian_Force_Densities(ds,Nb,xLag,yLag,springs,Lx,Ly);
    
else
    fx_springs = 0;%zeros(length(xLag),1); % No x-forces coming from springs
    fy_springs = fx_springs;               % No y-forces coming from springs
end



%----------------------------------------------------------------------
% Compute MASS PT FORCE DENSITIES (if there are mass points!)
%----------------------------------------------------------------------
if ( mass_Yes == 1)
    
    % Compute the Lagrangian MASSIVE PT force densities!
    [fx_mass, fy_mass, F_Mass] = give_Me_Mass_Lagrangian_Force_Densities(ds,xLag,yLag,masses); 
else
    fx_mass = 0;%zeros(length(xLag),1);; % No x-forces coming from mass points
    fy_mass = fx_mass;                   % No y-forces coming from mass points
    F_Mass = 0;                          % Dummy to pass along  
end


%--------------------------------------------------------------------
% Compute TARGET FORCE DENSITIES (if there are target points!)
%--------------------------------------------------------------------
if ( target_pts_Yes == 1)
    
    % Compute the Lagrangian TARGET force densities!
    [fx_target, fy_target] = give_Me_Target_Lagrangian_Force_Densities(ds,xLag,yLag,targets,Lx,Ly); 
    
else
    fx_target = 0;%zeros(length(xLag),1);; % No x-forces coming from target points
    fy_target = fx_target;                 % No y-forces coming from target points
end



%------------------------------------------------------------------------
% Compute BEAM (TORSIONAL SPRINGS) FORCE DENSITIES (if there are beams!)
%------------------------------------------------------------------------
if ( beams_Yes == 1 )

    % Compute the Lagrangian BEAM (TORSIONAL SPRINGS) force densities!
    [fx_beams, fy_beams] = give_Me_Beam_Lagrangian_Force_Densities(ds,Nb,xLag,yLag,beams,Lx,Ly);
    
else
    fx_beams = 0;%zeros(length(xLag),1); % No x-forces coming from beams
    fy_beams = fx_beams;                 % No y-forces coming from beams
end



%---------------------------------------------------------------------------------
% Compute BEAM (NON-INVARIANT) FORCE DENSITIES (if there are non-invariant beams!)
%---------------------------------------------------------------------------------
if ( nonInv_beams_Yes == 1 )

    % Compute the Lagrangian NON-INVARIANT BEAM force densities!
    [fx_nonInv_beams, fy_nonInv_beams] = give_Me_nonInv_Beam_Lagrangian_Force_Densities(ds,Nb,xLag,yLag,nonInv_beams,Lx,Ly);
    
else
    fx_nonInv_beams = 0;%zeros(length(xLag),1); % No x-forces coming from beams
    fy_nonInv_beams = fx_nonInv_beams;          % No y-forces coming from beams
end



%----------------------------------------------------------------------
% Compute DAMPED SPRING FORCE DENSITIES (if there are damped springs!)
%----------------------------------------------------------------------
if ( d_Springs_Yes == 1 )

    % Compute the Lagrangian DAMPED SPRING force densities!
    [fx_dSprings, fy_dSprings] = give_Me_Damped_Springs_Lagrangian_Force_Densities(ds,Nb,xLag,yLag,d_Springs,xLag_P,yLag_P,dt,Lx,Ly);
    
else
    fx_dSprings = 0;%zeros(length(xLag),1);;    %No x-forces coming from damped springs
    fy_dSprings = fx_dSprings;    %No y-forces coming from damped springs
end


%---------------------------------------------------------------------------------
% Compute GENERAL USER-DEFINED FORCE DENSITIES (if there is a user-defined force!)
%---------------------------------------------------------------------------------
if ( gen_force_Yes == 1 )

    % Compute the Lagrangian GENERAL USER force densities!
    [fx_genForce, fy_genForce] = give_Me_General_User_Defined_Force_Densities(ds,Nb,xLag,yLag,xLag_P,yLag_P,dt,current_time,general_force);
    
else
    fx_genForce = 0;%zeros(length(xLag),1);;  % No x-forces coming from general force model
    fy_genForce = fx_genForce;                % No y-forces coming from general force model
end

%--------------------------------------------------------------------------
% Compute COAGULATION MODEL AND FORCE DENSITIES (if there is coagulation!)
%--------------------------------------------------------------------------
if ( coagulation_Yes == 1 )

    % Compute the Lagrangian COAGULATION force densities!
    [fx_coag, fy_coag, aggregate_list] = give_Me_Coagulation_Force_Densities(Nb,xLag,yLag,coagulation,aggregate_list,Lx,Ly);
    
else
    fx_coag = 0;%zeros(length(xLag),1);  % No x-forces coming from coagulation
    fy_coag = fx_coag;                   % No y-forces coming from coagulation
end



%-----------------------------------------------------
% SUM TOTAL FORCE DENSITY! %
%-----------------------------------------------------
fx = fx_springs + fx_target + fx_beams + fx_nonInv_beams + fx_muscles + fx_muscles3 + fx_mass + fx_dSprings + fx_genForce + fx_coag;
fy = fy_springs + fy_target + fy_beams + fy_nonInv_beams + fy_muscles + fy_muscles3 + fy_mass + fy_dSprings + fy_genForce + fy_coag;


%-----------------------------------------------------
% Save Poro-Elastic Forces, if poroelastic elements %
%-----------------------------------------------------
if poroelastic_Yes  
    F_Poro(:,1) = fx_springs(poroelastic_info(:,1)) + fx_nonInv_beams(poroelastic_info(:,1));
    F_Poro(:,2) = fy_springs(poroelastic_info(:,1)) + fy_nonInv_beams(poroelastic_info(:,1));
else
    F_Poro = 0;
end


%----------------------------------
% SAVE LAGRANGIAN FORCES
%----------------------------------
F_Lag(:,1) = fx;
F_Lag(:,2) = fy;
    

%-----------------------------------------
% Give me delta-function approximations!
%-----------------------------------------
[delta_X, delta_Y] = give_Me_Delta_Function_Approximations_For_Force_Calc(x,y,grid_Info,xLag,yLag);


%------------------------------------------------------------
% Transform the force density vectors into diagonal matrices
%------------------------------------------------------------
%
if flag_Geo_Connect % Compute -actual- distances between 
                    % neighboring Lagrangian pts in geometry
    
    for i=1:Nb
        
        %-----------------------------------------------------------
        % Compute real distances btwn LAG_i and Attached Points
        %         and sum distances together
        %-----------------------------------------------------------
        
        % find all Lag Pts that are Lag_i's geometric neighbors
        indsVec = find(geo_Connect_MAT(:,1)==i);
    
        % loop to (1) distances between Lag_i and each neighbor
        %         (2) sum all such distances together
        ds_sum = 0;
        for j=1:length(indsVec)
            id1 = geo_Connect_MAT( indsVec(j),1);
            id2 = geo_Connect_MAT( indsVec(j),2);
            ds_sum = ds_sum + sqrt( ( xLag(id1,1)-xLag(id2,1) )^2 + ( yLag(id1,1)-yLag(id2,1) )^2 );     
        end
        
        %---------------------------------------------------------------
        %   ORIGINAL BEFORE MORE EFFICIENT MATRIX MULTIPLICATION BELOW!
        % Note: (1) 0.5 coming from Trapezoid Rule
        %       (2) ds_sum is sum of distances to neighboring pts
        %---------------------------------------------------------------
        fx(i,1) = 0.5*fx(i,1)*ds_sum; 
        fy(i,1) = 0.5*fy(i,1)*ds_sum; 
        
    end
    
else
    
    %-----------------------------------------
    % Use Peskin's Constant 'ds' assumption
    %-----------------------------------------
    %for i=1:Nb
    %   fxds(i,i) = fx(i,1)*ds; 
    %   fyds(i,i) = fy(i,1)*ds;
    %end
    
    %-----------------------------------------
    % Uses Peskin's Constant 'ds' assumption
    % (vectorized instead of for-loop above)
    %-----------------------------------------
    %fxds = diag(fx(:,1)*ds);
    %fyds = diag(fy(:,1)*ds);
    
    %---------------------------------------------------
    %   Uses Peskin's Constant 'ds' assumption
    % (Using efficient diag-matrix multiplier below)
    %---------------------------------------------------
    fx = fx(:,1)*ds;
    fy = fy(:,1)*ds;
    
end


%-----------------------------------------------------------------------
%   ORIGINAL --> DOESN'T USE EFFICIENT MATRIX MULTIPLICATION W/ 
%                INVOLVING DIAGONAL MATRIX
%
% Find Eulerian forces on grids by approximating the line integral, 
%       F(x,y) = int{ f(s) delta(x - xLag(s)) delta(y - yLag(s)) ds }
%-----------------------------------------------------------------------
% fxds = diag(fx(:,1)*ds);
% fyds = diag(fy(:,1)*ds);
% Fx = delta_Y * fxds * delta_X;
% Fy = delta_Y * fyds * delta_X;


%-----------------------------------------------------------------------
%    USES MORE EFFICIENT MATRIX MULTIPLICATION W/ A DIAGONAL MATRIX
%
% Find Eulerian forces on grids by approximating the line integral, 
%       F(x,y) = int{ f(s) delta(x - xLag(s)) delta(y - yLag(s)) ds }
%            eg,
%                 A*diag(vec) = ( (diag(vec)).*A' )'
%             (note the transposes in the above line)
%-----------------------------------------------------------------------
Fx = ( ( fx(:,1) ).*delta_Y' )' * delta_X;
Fy = ( ( fy(:,1) ).*delta_Y' )' * delta_X;










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION computes the Lagrangian SPRING Force Densities.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fx, fy] = give_Me_Spring_Lagrangian_Force_Densities(ds,Nb,xLag,yLag,springs,Lx,Ly)


Nsprings = length(springs(:,1));  % # of Springs
sp_1 = springs(:,1);              % Initialize storage for LEADER NODE Spring Connection
sp_2 = springs(:,2);              % Initialize storage for FOLLOWER NODE Spring Connection
K_Vec = springs(:,3);             % Stores spring stiffness associated with each spring
RL_Vec = springs(:,4);            % Stores spring resting length associated with each spring
alpha_pow = springs(:,5);         % Degree of linearity (1=linear, >1 = non-linear)

fx = zeros(Nb,1);                 % Initialize storage for x-forces
fy = fx;                          % Initialize storage for y-forces


%-----------------------------------------------
% LOOP THROUGH ALL LAGRANGIAN SPRINGS
%-----------------------------------------------
for i=1:Nsprings
    
    id_Leader = sp_1(i);          % Leader Node index
    id_Follower = sp_2(i);        % Follower Node index
    k_Spring = K_Vec(i);          % Spring stiffness of i-th spring
    L_r = RL_Vec(i);              % Resting length of i-th spring
    alpha = alpha_pow(i);         % Degree of linearity of i-th spring
 
    
    dx = xLag(id_Follower) - xLag(id_Leader); % x-Distance btwn follower and leader node
    dy = yLag(id_Follower) - yLag(id_Leader); % y-Distance btwn follower and leader node

    %
    % TESTING FOR LAG PT. PASSED THRU BNDRY; MAY NEED TO CHANGE TOLERANCE HERE, DEPENDENT ON APPLICATION
    %
    if abs(dx) > Lx/2
        dx = sign(dx)*( Lx - sign(dx)*dx );
    end
    
    if abs(dy) > Ly/2
        dy = sign(dy)*( Ly - sign(dy)*dy );
    end
    
    sF_x = 0.5*(alpha+1) * k_Spring * ( sqrt( dx^2 + dy^2 ) - L_r )^(alpha) * ( dx / sqrt(dx^2+dy^2) );
    sF_y = 0.5*(alpha+1) * k_Spring * ( sqrt( dx^2 + dy^2 ) - L_r )^(alpha) * ( dy / sqrt(dx^2+dy^2) );
    
    fx(id_Leader,1) = fx(id_Leader,1) + sF_x;  % Sum total forces for node, i in x-direction (this is LEADER node for this spring)
    fy(id_Leader,1) = fy(id_Leader,1) + sF_y;  % Sum total forces for node, i in y-direction (this is LEADER node for this spring)
    
    fx(id_Follower,1) = fx(id_Follower,1) - sF_x;  % Sum total forces for node, i in x-direction (this is FOLLOWER node for this spring)
    fy(id_Follower,1) = fy(id_Follower,1) - sF_y;  % Sum total forces for node, i in y-direction (this is FOLLOWER node for this spring)

    
end

% MIGHT NOT NEED THESE!
%fx = fx/ds^2;
%fy = fy/ds^2;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION computes the Lagrangian **DAMPED** SPRING Force Densities .
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fx, fy] = give_Me_Damped_Springs_Lagrangian_Force_Densities(ds,Nb,xLag,yLag,d_Springs,xLag_P,yLag_P,dt,Lx,Ly)


Ndsprings = length(d_Springs(:,1));  % # of DAMPED Springs
sp_1 = d_Springs(:,1);               % Initialize storage for LEADER NODE Spring Connection
sp_2 = d_Springs(:,2);               % Initialize storage for FOLLOWER NODE Spring Connection
K_Vec = d_Springs(:,3);              % Stores spring stiffness associated with each spring
RL_Vec = d_Springs(:,4);             % Stores spring resting length associated with each spring
b_Vec = d_Springs(:,5);              % Damping coefficient

fx = zeros(Nb,1);                    % Initialize storage for x-forces
fy = fx;                             % Initialize storage for y-forces


%-----------------------------------------------
% LOOP THROUGH ALL DAMPED LAGRANGIAN SPRINGS
%-----------------------------------------------
for i=1:Ndsprings
    
    id_Leader = sp_1(i);          % Leader Node index
    id_Follower = sp_2(i);        % Follower Node index
    k_Spring = K_Vec(i);          % Spring stiffness of i-th spring
    L_r = RL_Vec(i);              % Resting length of i-th spring
    b = b_Vec(i);                 % Damping Coefficient
    
    dx = xLag(id_Follower) - xLag(id_Leader);      % x-Distance btwn follower and leader node
    dy = yLag(id_Follower) - yLag(id_Leader);      % y-Distance btwn follower and leader node
    
    %
    % TESTING FOR LAG PT. PASSED THRU BNDRY; MAY NEED TO CHANGE TOLERANCE HERE, DEPENDENT ON APPLICATION
    %
    if abs(dx) > Lx/2
       dx = sign(dx)*( Lx - sign(dx)*dx );
    end
    
    if abs(dy) > Lx/2
        dy = sign(dy)*( Ly - sign(dy)*dy );
    end
    
    dV_x = ( xLag(id_Leader) - xLag_P(id_Leader) ); % dt*(x-Velocity) between current and prev. steps
    dV_y = ( yLag(id_Leader) - yLag_P(id_Leader) ); % dt*(y-Velocity) between current and prev. steps
    
    %
    % TESTING FOR LAG PT. PASSED THRU BNDRY; MAY NEED TO CHANGE TOLERANCE HERE, DEPENDENT ON APPLICATION
    %
    if abs(dV_x) > Lx/2
       dV_x = sign(dV_x)*( Lx - sign(dV_x)*dV_x );
    end
    
    if abs(dV_y) > Lx/2
        dV_y = sign(dV_y)*( Ly - sign(dV_y)*dV_y );
    end
    
    
    dV_x = ( dV_x ) / dt ; % x-Velocity between current and prev. steps
    dV_y = ( dV_y ) / dt ; % y-Velocity between current and prev. steps
    
    
    sF_x = k_Spring * ( sqrt( dx^2 + dy^2 ) - L_r ) * ( dx / sqrt(dx^2+dy^2) ) - b*dV_x; %added negative for testing
    sF_y = k_Spring * ( sqrt( dx^2 + dy^2 ) - L_r ) * ( dy / sqrt(dx^2+dy^2) ) - b*dV_y;
    
    fx(id_Leader,1) = fx(id_Leader,1) + sF_x ;  % Sum total forces for node, i in x-direction (this is LEADER node for this spring)
    fy(id_Leader,1) = fy(id_Leader,1) + sF_y ;  % Sum total forces for node, i in y-direction (this is LEADER node for this spring)
    
    fx(id_Follower,1) = fx(id_Follower,1) - sF_x ; % Sum total forces for node, i in x-direction (this is FOLLOWER node for this spring)
    fy(id_Follower,1) = fy(id_Follower,1) - sF_y ; % Sum total forces for node, i in y-direction (this is FOLLOWER node for this spring)

    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION computes the Lagrangian MUSCLE Force Densities for LENGTH-TENSION/FORCE-VELOCITY MUSCLE MODEL.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fx,fy] = give_ElectroPhys_Ca_Muscle_Force_Densities(Nb,xLag,yLag,xLag_P,yLag_P,muscles,current_time,dt,electro_potential)


Nmuscles = length(muscles(:,1));  % # of Muscles
m_1 = muscles(:,1);               % Initialize storage for LEADER NODE Muscle Connection
m_2 = muscles(:,2);               % Initialize storage for FOLLOWER NODE Muscle Connection
LFO_Vec = muscles(:,3);           % Stores length for max. muscle tension
SK_Vec = muscles(:,4);            % Stores muscle constant
a_Vec = muscles(:,5);             % Stores Hill Parameter, a
b_Vec = muscles(:,6);             % Stores Hill Parameter, b
FMAX_Vec = muscles(:,7);          % Stores Force-Maximum for Muscle

fx = zeros(Nb,1);                 % Initialize storage for x-forces
fy = fx;                          % Initialize storage for y-forces

ct = current_time/dt;             % gives count of time-steps for simulation

%-----------------------------------------------
% LOOP THROUGH ALL LAGRANGIAN MUSCLES
%-----------------------------------------------
for i=1:Nmuscles
    
    id_Leader = m_1(i);          % Leader Node index for i-th muscle
    id_Follower = m_2(i);        % Follower Node index for i-th muscle
    LFO = LFO_Vec(i);            % Length for max. muscle tension for i-th muscle
    sk = SK_Vec(i);              % Muscle constant for i-th muscle
    a = a_Vec(i);                % Hill parameter, a, for i-th muscle
    b = b_Vec(i);                % Hill parameter, b, for i-th muscle
    Fmax = FMAX_Vec(i);          % Force-Maximum for i-th muscle
    
    %xPt = xLag( id_Leader );     % x-Pt of interest at the moment to drive muscle contraction
    
    dx = xLag(id_Follower) - xLag(id_Leader); % x-Distance btwn follower and leader node
    dy = yLag(id_Follower) - yLag(id_Leader); % y-Distance btwn follower and leader node
    LF = sqrt( dx^2 + dy^2 );                 % Euclidean DISTANCE between leader and follower node
    
    
    dx_P = xLag_P(id_Follower) - xLag_P(id_Leader); % x-Distance btwn follower and leader node
    dy_P = yLag_P(id_Follower) - yLag_P(id_Leader); % y-Distance btwn follower and leader node
    LF_P = sqrt( dx_P^2 + dy_P^2 );                 % Euclidean DISTANCE between leader and follower node
    
    v =  abs(LF-LF_P)/dt;        % How fast the muscle is contracting/expanding 

    % Find actual muscle activation magnitude
    Fm = give_Ca_ElectroPhys_Muscle_Activation(v,LF,LFO,sk,a,b,Fmax,ct,id_Leader,electro_potential);
    
    mF_x = Fm*(dx/LF);           % cos(theta) = dx / LF;
    mF_y = Fm*(dy/LF);           % sin(theta) = dy / LF;
    
    fx(id_Leader,1) = fx(id_Leader,1) + mF_x;  % Sum total forces for node, i in x-direction (this is LEADER node for this spring)
    fy(id_Leader,1) = fy(id_Leader,1) + mF_y;  % Sum total forces for node, i in y-direction (this is LEADER node for this spring)
    
    fx(id_Follower,1) = fx(id_Follower,1) - mF_x;  % Sum total forces for node, i in x-direction (this is FOLLOWER node for this spring)
    fy(id_Follower,1) = fy(id_Follower,1) - mF_y;  % Sum total forces for node, i in y-direction (this is FOLLOWER node for this spring)

    
end











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION computes the Lagrangian MUSCLE Force Densities for LENGTH-TENSION/FORCE-VELOCITY MUSCLE MODEL.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fx,fy] = give_Muscle_Force_Densities(Nb,xLag,yLag,xLag_P,yLag_P,muscles,current_time,dt)


Nmuscles = length(muscles(:,1));  % # of Muscles
m_1 = muscles(:,1);               % Initialize storage for LEADER NODE Muscle Connection
m_2 = muscles(:,2);               % Initialize storage for FOLLOWER NODE Muscle Connection
LFO_Vec = muscles(:,3);           % Stores length for max. muscle tension
SK_Vec = muscles(:,4);            % Stores muscle constant
a_Vec = muscles(:,5);             % Stores Hill Parameter, a
b_Vec = muscles(:,6);             % Stores Hill Parameter, b
FMAX_Vec = muscles(:,7);          % Stores Force-Maximum for Muscle

fx = zeros(Nb,1);                 % Initialize storage for x-forces
fy = fx;                          % Initialize storage for y-forces

%-----------------------------------------------
% LOOP THROUGH ALL LAGRANGIAN MUSCLES
%-----------------------------------------------
for i=1:Nmuscles
    
    id_Leader = m_1(i);          % Leader Node index for i-th muscle
    id_Follower = m_2(i);           % Follower Node index for i-th muscle
    LFO = LFO_Vec(i);            % Length for max. muscle tension for i-th muscle
    sk = SK_Vec(i);              % Muscle constant for i-th muscle
    a = a_Vec(i);                % Hill parameter, a, for i-th muscle
    b = b_Vec(i);                % Hill parameter, b, for i-th muscle
    Fmax = FMAX_Vec(i);          % Force-Maximum for i-th muscle
    
    xPt = xLag( id_Leader );     % x-Pt of interest at the moment to drive muscle contraction
    
    dx = xLag(id_Follower) - xLag(id_Leader); % x-Distance btwn follower and leader node
    dy = yLag(id_Follower) - yLag(id_Leader); % y-Distance btwn follower and leader node
    LF = sqrt( dx^2 + dy^2 );              % Euclidean DISTANCE between leader and follower node
    
    
    dx_P = xLag_P(id_Follower) - xLag_P(id_Leader); % x-Distance btwn follower and leader node
    dy_P = yLag_P(id_Follower) - yLag_P(id_Leader); % y-Distance btwn follower and leader node
    LF_P = sqrt( dx_P^2 + dy_P^2 );              % Euclidean DISTANCE between leader and follower node
    
    v =  abs(LF-LF_P)/dt;        % How fast the muscle is contracting/expanding 

    % Find actual muscle activation magnitude
    Fm = give_Muscle_Activation(v,LF,LFO,sk,a,b,Fmax,current_time,xPt,xLag);
    
    mF_x = Fm*(dx/LF);           % cos(theta) = dx / LF;
    mF_y = Fm*(dy/LF);           % sin(theta) = dy / LF;
    
    fx(id_Leader,1) = fx(id_Leader,1) + mF_x;  % Sum total forces for node, i in x-direction (this is LEADER node for this spring)
    fy(id_Leader,1) = fy(id_Leader,1) + mF_y;  % Sum total forces for node, i in y-direction (this is LEADER node for this spring)
    
    fx(id_Follower,1) = fx(id_Follower,1) - mF_x;    % Sum total forces for node, i in x-direction (this is FOLLOWER node for this spring)
    fy(id_Follower,1) = fy(id_Follower,1) - mF_y;    % Sum total forces for node, i in y-direction (this is FOLLOWER node for this spring)

    
end












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes the Lagrangian MUSCLE Force Densities for 3-ELEMENT
%           HILL MODEL! ( coupled w/ force-velocity/length-tension model for
%           contractile element!)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fx,fy] = give_3_Element_Muscle_Force_Densities(Nb,xLag,yLag,xLag_P,yLag_P,muscles3,current_time,dt)


Nmuscles = length(muscles3(:,1));  % # of Muscles
m_1 = muscles3(:,1);               % Initialize storage for LEADER NODE Muscle Connection
m_2 = muscles3(:,2);               % Initialize storage for FOLLOWER NODE Muscle Connection
LFO_Vec = muscles3(:,3);           % Stores length for max. muscle tension
SK_Vec = muscles3(:,4);            % Stores muscle constant
a_Vec = muscles3(:,5);             % Stores Hill Parameter, a
b_Vec = muscles3(:,6);             % Stores Hill Parameter, b
FMAX_Vec = muscles3(:,7);          % Stores Force-Maximum for Muscle
kSpr_Vec = muscles3(:,8);          % Stores Spring Coeffs
alpha_pow= muscles3(:,9);          % Stores deg. of non-linearity for springs

fx = zeros(Nb,1);                 % Initialize storage for x-forces
fy = fx;                          % Initialize storage for y-forces

%-----------------------------------------------
% LOOP THROUGH ALL LAGRANGIAN MUSCLES
%-----------------------------------------------
for i=1:Nmuscles
    
    id_Leader = m_1(i);          % Leader Node index for i-th muscle
    id_Follower = m_2(i);        % Follower Node index for i-th muscle
    LFO = LFO_Vec(i);            % Length for max. muscle tension for i-th muscle
    sk = SK_Vec(i);              % Muscle constant for i-th muscle
    a = a_Vec(i);                % Hill parameter, a, for i-th muscle
    b = b_Vec(i);                % Hill parameter, b, for i-th muscle
    Fmax = FMAX_Vec(i);          % Force-Maximum for i-th muscle
    kSpr = kSpr_Vec(i);          % Spring coefficient for PARALLEL ELEMENT NL-spring
    alpha = alpha_pow(i);        % Degree of linearity of PARALLEL ELEMENT i-th spring

    
    xPt = xLag( id_Leader );     % x-Pt of interest at the moment to drive muscle contraction
    
    dx = xLag(id_Follower) - xLag(id_Leader); % x-Distance btwn follower and leader node
    dy = yLag(id_Follower) - yLag(id_Leader); % y-Distance btwn follower and leader node
    LF = sqrt( dx^2 + dy^2 );                 % Euclidean DISTANCE between leader and follower node
    
    
    dx_P = xLag_P(id_Follower) - xLag_P(id_Leader); % x-Distance btwn follower and leader node
    dy_P = yLag_P(id_Follower) - yLag_P(id_Leader); % y-Distance btwn follower and leader node
    LF_P = sqrt( dx_P^2 + dy_P^2 );                 % Euclidean DISTANCE between leader and follower node
    
    v =  abs(LF-LF_P)/dt;        % How fast the muscle is contracting/expanding 
    
    % Find actual muscle activation magnitude for CONTRACTILE ELEMENT
    [Fm,on_Parallel,PE_Coeff,af_Val] = give_3_Element_Muscle_Activation(v,LF,LFO,sk,a,b,Fmax,current_time,xPt,xLag,i);

    
    % Compute muscle force from CONTRACTILE ELEMENT in each direction
    mF_x = Fm*(dx/LF);           % cos(theta) = dx / LF;
    mF_y = Fm*(dy/LF);           % sin(theta) = dy / LF;
    
    
    % Find muscle force from SERIES ELEMENT in each direction
    bDamp = 1.0; 
    dV_x = ( xLag(id_Leader) - xLag_P(id_Leader) )/dt;      % Compute velocity gradient terms for damping
    dV_y = ( yLag(id_Leader) - yLag_P(id_Leader) )/dt;      % Compute velocity gradient terms for damping
    sF_SE_x = af_Val * kSpr * ( (LFO-LF) - LFO ) * ( dx / sqrt(dx^2+dy^2) ) + bDamp*dV_x; % Note: L_con (LF) + L_ser = L_tot = LFO
    sF_SE_y = af_Val * kSpr * ( (LFO-LF) - LFO ) * ( dy / sqrt(dx^2+dy^2) ) + bDamp*dV_y;
    
    
    
    if on_Parallel == 1
        % Find force from series element (tendons)
        sF_PE_x = PE_Coeff * 0.5*(alpha+1) * kSpr * ( sqrt( dx^2 + dy^2 ) - LFO )^(alpha) * ( dx / sqrt(dx^2+dy^2) );
        sF_PE_y = PE_Coeff * 0.5*(alpha+1) * kSpr * ( sqrt( dx^2 + dy^2 ) - LFO )^(alpha) * ( dy / sqrt(dx^2+dy^2) );
    else
        % Parallel Element is off
        sF_PE_x = 0;
        sF_PE_y = 0;
    end
    
    
    fx(id_Leader,1) = fx(id_Leader,1) + mF_x + sF_SE_x + sF_PE_x;  % Sum total forces for node, i in x-direction (this is LEADER node for this spring)
    fy(id_Leader,1) = fy(id_Leader,1) + mF_y + sF_SE_y + sF_PE_y;  % Sum total forces for node, i in y-direction (this is LEADER node for this spring)
    
    fx(id_Follower,1) = fx(id_Follower,1) - mF_x - sF_SE_x - sF_PE_x; % Sum total forces for node, i in x-direction (this is FOLLOWER node for this spring)
    fy(id_Follower,1) = fy(id_Follower,1) - mF_y - sF_SE_y - sF_PE_y; % Sum total forces for node, i in y-direction (this is FOLLOWER node for this spring)


end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION computes the Lagrangian BEAM (NON-INVARIANT) Force Densities 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fx, fy] = give_Me_nonInv_Beam_Lagrangian_Force_Densities(ds,Nb,xLag,yLag,beams,Lx,Ly)

Nbeams = length(beams(:,1));     % # of Beams
pts_1 = beams(:,1);              % Initialize storage for 1ST NODE for BEAM
pts_2 = beams(:,2);              % Initialize storage for MIDDLE NODE (2ND Node) for BEAM
pts_3 = beams(:,3);              % Initialize storage for 3RD NODE for BEAM
K_Vec = beams(:,4);              % Stores beam stiffness associated with each beam
CX_Vec = beams(:,5);             % Stores beam curvature in x-direction
CY_Vec = beams(:,6);             % Stores beam curvature in y-direction

fx = zeros(Nb,1);                % Initialize storage for x-forces
fy = fx;                         % Initialize storage for y-forces

%-----------------------------------------------
% LOOP THROUGH ALL LAGRANGIAN NON-INV BEAMS
%-----------------------------------------------
for i=1:Nbeams
    
    id_1 = pts_1(i);          % 1ST Node index
    id_2 = pts_2(i);          % (MIDDLE) 2nd Node index -> index that gets force applied to it!
    id_3 = pts_3(i);          % 3RD Node index
    k_Beam = K_Vec(i);        % Beam stiffness of i-th spring
    Cx = CX_Vec(i) ;          % x-Curvature of the beam between these three nodes 
    Cy = CY_Vec(i) ;          % y-Curvature of the beam between these three nodes 

    Xp = xLag(id_1);          % xPt of 1ST Node Pt. in beam
    Xq = xLag(id_2);          % xPt of 2ND (MIDDLE) Node Pt. in beam
    Xr = xLag(id_3);          % xPt of 3RD Node Pt. in beam
    
    Yp = yLag(id_1);          % yPt of 1ST Node Pt. in beam
    Yq = yLag(id_2);          % yPt of 2ND (MIDDLE) Node Pt. in beam
    Yr = yLag(id_3);          % yPt of 3RD Node Pt. in beam
    
    % Checks if Lag. Pts. have passed through the boundary and translates appropriately
    [Xp,Xq,Xr] = check_If_Beam_Points_Pass_Through_Boundary(ds,Lx,Xp,Xq,Xr);
    [Yp,Yq,Yr] = check_If_Beam_Points_Pass_Through_Boundary(ds,Ly,Yp,Yq,Yr);
    
    % CALCULATE BENDING IN X
    fx(id_3,1) = fx(id_3,1) -   k_Beam*( Xr - 2*Xq + Xp - Cx);
    fx(id_1,1) = fx(id_1,1) -   k_Beam*( Xr - 2*Xq + Xp - Cx);
    fx(id_2,1) = fx(id_2,1) + 2*k_Beam*( Xr - 2*Xq + Xp - Cx );
    
    % CALCULATE BENDING IN Y
    fy(id_3,1) = fy(id_3,1) -   k_Beam*( Yr - 2*Yq + Yp - Cy );
    fy(id_1,1) = fy(id_1,1) -   k_Beam*( Yr - 2*Yq + Yp - Cy );
    fy(id_2,1) = fy(id_2,1) + 2*k_Beam*( Yr - 2*Yq + Yp - Cy );
    
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION computes the Lagrangian BEAM (TORSIONAL SPRING) Force Densities 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fx, fy] = give_Me_Beam_Lagrangian_Force_Densities(ds,Nb,xLag,yLag,beams,Lx,Ly)


Nbeams = length(beams(:,1));     % # of Beams
pts_1 = beams(:,1);              % Initialize storage for 1ST NODE for BEAM
pts_2 = beams(:,2);              % Initialize storage for MIDDLE NODE (2ND Node) for BEAM
pts_3 = beams(:,3);              % Initialize storage for 3RD NODE for BEAM
K_Vec = beams(:,4);              % Stores spring stiffness associated with each spring
C_Vec = beams(:,5);              % Stores spring resting length associated with each spring

fx = zeros(Nb,1);                % Initialize storage for x-forces
fy = fx;                         % Initialize storage for y-forces

%--------------------------------------------------------
% LOOP THROUGH ALL LAGRANGIAN BEAMS (Torsional Springs)
%--------------------------------------------------------
for i=1:Nbeams
    
    id_1 = pts_1(i);          % 1ST Node index
    id_2 = pts_2(i);          % (MIDDLE) 2nd Node index -> index that gets force applied to it!
    id_3 = pts_3(i);          % 3RD Node index
    k_Beam = K_Vec(i);        % Beam stiffness of i-th spring
    C = C_Vec(i) ;            % Curvature of the beam between these three nodes 
    
    Xp = xLag(id_1);          % xPt of 1ST Node Pt. in beam
    Xq = xLag(id_2);          % xPt of 2ND (MIDDLE) Node Pt. in beam
    Xr = xLag(id_3);          % xPt of 3RD Node Pt. in beam
    
    Yp = yLag(id_1);          % yPt of 1ST Node Pt. in beam
    Yq = yLag(id_2);          % yPt of 2ND (MIDDLE) Node Pt. in beam
    Yr = yLag(id_3);          % yPt of 3RD Node Pt. in beam
    
    % Checks if Lag. Pts. have passed through the boundary and translates appropriately
    [Xp,Xq,Xr] = check_If_Beam_Points_Pass_Through_Boundary(ds,Lx,Xp,Xq,Xr);
    [Yp,Yq,Yr] = check_If_Beam_Points_Pass_Through_Boundary(ds,Ly,Yp,Yq,Yr);
   
    % Compute Cross-Product
    cross_prod = (Xr-Xq)*(Yq-Yp) - (Yr-Yq)*(Xq-Xp);
    
    % FORCES FOR LEFT NODE
    bF_x_L = -k_Beam * ( cross_prod - C ) * ( Yr-Yq );
    bF_y_L =  k_Beam * ( cross_prod - C ) * ( Xr-Xq );
    
    % FORCES FOR MIDDLE NODE
    bF_x_M =  k_Beam * ( cross_prod - C ) * (  (Yq-Yp) + (Yr-Yq) );
    bF_y_M = -k_Beam * ( cross_prod - C ) * (  (Xr-Xq) + (Xq-Xp) );
    
    % FORCES FOR RIGHT NODE
    bF_x_R = -k_Beam * ( cross_prod - C ) * ( Yq-Yp );
    bF_y_R =  k_Beam * ( cross_prod - C ) * ( Xq-Xp );
    
    fx(id_1,1) = fx(id_1,1) - bF_x_L;  % Sum total forces for left node, in x-direction (this is LEFT node for this beam)
    fy(id_1,1) = fy(id_1,1) - bF_y_L;  % Sum total forces for left node, in y-direction (this is LEFT node for this beam)
    
    fx(id_2,1) = fx(id_2,1) + bF_x_M;  % Sum total forces for middle node, in x-direction (this is MIDDLE node for this beam)
    fy(id_2,1) = fy(id_2,1) + bF_y_M;  % Sum total forces for middle node, in y-direction (this is MIDDLE node for this beam)
    
    fx(id_3,1) = fx(id_3,1) - bF_x_R;  % Sum total forces for right node, in x-direction (this is RIGHT node for this beam)
    fy(id_3,1) = fy(id_3,1) - bF_y_R;  % Sum total forces for right node, in y-direction (this is RIGHT node for this beam)
    
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes the MASS Force Densities! 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fx_mass, fy_mass, F_Mass] = give_Me_Mass_Lagrangian_Force_Densities(ds,xLag,yLag,masses)

IDs = masses(:,1);                 % Stores Lag-Pt IDs in col vector
xPts= masses(:,2);                 % Original x-Values of x-Mass Pts.
yPts= masses(:,3);                 % Original y-Values of y-Mass Pts.
kStiffs = masses(:,4);             % Stores "spring" stiffness parameter

N_masses = length(IDs);            % # of target points!

fx = zeros(length(xLag),1);         % Initialize storage for x-force density from TARGET PTS
fy = fx;                            % Initialize storage for y-force density from TARGET PTS

%-------------------------------------------
% LOOP THROUGH ALL LAGRANGIAN MASS POINTS 
%-------------------------------------------
for i=1:N_masses
   
    fx(IDs(i),1) = fx(IDs(i),1) + kStiffs(i)*( xPts(i) - xLag(IDs(i)) );
    fy(IDs(i),1) = fy(IDs(i),1) + kStiffs(i)*( yPts(i) - yLag(IDs(i)) ); 
   
end

fx_mass = fx;
fy_mass = fy;

F_Mass(:,1) = fx;  % Store for updating massive boundary pts
F_Mass(:,2) = fy;  % Store for updating massive boundary pts

% MIGHT NOT NEED THESE!
%fx_mass = fx/ds^2;
%fy_mass = fy/ds^2;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes the Target-Pt Force Densities! 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fx_target, fy_target] = give_Me_Target_Lagrangian_Force_Densities(ds,xLag,yLag,targets,Lx,Ly)

IDs = targets(:,1);                 % Stores Lag-Pt IDs in col vector
xPts= targets(:,2);                 % Original x-Values of x-Target Pts.
yPts= targets(:,3);                 % Original y-Values of y-Target Pts.
kStiffs = targets(:,4);             % Stores Target Stiffnesses 

N_targets = length( targets(:,1) ); % # of target points!

fx = zeros(length(xLag),1);         % Initialize storage for x-force density from TARGET PTS
fy = fx;                            % Initialize storage for y-force density from TARGET PTS


%----------------------------------------------------------
% Original...same CPU time roughly as...
%      --> NOT-initializing other data structures
%      --> Vectorizing the computation
%----------------------------------------------------------
for i=1:N_targets
   
    dx = xPts(i) - xLag(IDs(i)); % x-Distance btwn Lag Pt. and Virtual pt
    dy = yPts(i) - yLag(IDs(i)); % y-Distance btwn Lag Pt. and Virtual pt
    
    %
    % TESTING FOR LAG PT. PASSED THRU BNDRY; MAY NEED TO CHANGE TOLERANCE HERE, DEPENDENT ON APPLICATION
    %
    if abs(dx) > Lx/2
        dx = sign(dx)*( Lx - sign(dx)*dx );
    end
    
    if abs(dy) > Ly/2
        dy = sign(dy)*( Ly - sign(dy)*dy );
    end  
    
    fx(IDs(i),1) = fx(IDs(i),1) + kStiffs(i)*( dx );
    fy(IDs(i),1) = fy(IDs(i),1) + kStiffs(i)*( dy ); 
   
end

%---------------------------------------------------------------------
% Original...but without defining other data-structures, as above
%               --> NO REAL CPU TIME SAVINGS...
%---------------------------------------------------------------------
% for i=1:N_targets
%    
%     dx = targets(i,2) - xLag(targets(i,1)); % x-Distance btwn Lag Pt. and Virtual pt
%     dy = targets(i,3) - yLag(targets(i,1)); % y-Distance btwn Lag Pt. and Virtual pt
%     
%     %
%     % TESTING FOR LAG PT. PASSED THRU BNDRY; MAY NEED TO CHANGE TOLERANCE HERE, DEPENDENT ON APPLICATION
%     %
%     if abs(dx) > Lx/2
%         dx = sign(dx)*( Lx - sign(dx)*dx );
%     end
%     
%     if abs(dy) > Ly/2
%         dy = sign(dy)*( Ly - sign(dy)*dy );
%     end  
%     
%     fx(targets(i,1),1) = fx(targets(i,1),1) + targets(i,4)*( dx );
%     fy(targets(i,1),1) = fy(targets(i,1),1) + targets(i,4)*( dy ); 
%    
% end

%----------------------------------------------------------
% Vectorized Version...same as non-vectorized version 
%               --> NO REAL CPU TIME SAVINGS...
%----------------------------------------------------------
% fx(IDs,1) = fx(IDs,1) + kStiffs .* ( xPts - xLag(IDs) );
% fy(IDs,1) = fy(IDs,1) + kStiffs .* ( yPts - yLag(IDs) );


fx_target = fx;
fy_target = fy;

% MIGHT NOT NEED THESE!
%fx_target = fx/ds^2;
%fy_target = fy/ds^2;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION computes the Delta-Function Approximations 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [delta_X, delta_Y] = give_Me_Delta_Function_Approximations_For_Force_Calc(x,y,grid_Info,xLag,yLag)

% Grid Info
Nx =   grid_Info(1);
Ny =   grid_Info(2);
Lx =   grid_Info(3);
Ly =   grid_Info(4);
dx =   grid_Info(5);
dy =   grid_Info(6);
supp = grid_Info(7);
Nb =   grid_Info(8);


%--------------------------------------------------------------
% Find the indices of the points (xi, yj) where the 
%           1D delta functions are non-zero in EULERIAN FRAME
%--------------------------------------------------------------
indX = give_1D_NonZero_Delta_Indices(xLag, Nx, dx, supp);
indY = give_1D_NonZero_Delta_Indices(yLag, Ny, dy, supp)';


%--------------------------------------------------------------
% Matrix of possible indices, augmented by "supp"-copies to 
%           perform subtractions later in LAGRANGIAN FRAME
%--------------------------------------------------------------
indLagAux = (1:1:Nb)';


%----------------------------------------------------
% ORIGINAL -> doesn't seem slower than repmat method
%----------------------------------------------------
ind_Lag = zeros(Nb,supp);
for i=1:supp
   ind_Lag(:,i)=indLagAux;
end

%--------------------------------------------------------------
% TRY USING REPMAT TO GET AROUND INITIALIZATION ISSUE ABOVE!
%       --> Doesn't make a difference
%       --> Tried case w/ 4096 Lag. Pts...
%--------------------------------------------------------------
%ind_Lag = repmat(indLagAux,1,supp);


%--------------------------------------------------------------
% Compute distance btwn Eulerian Pts and Lagrangian Pts 
%                       by passing correct indices for each
%--------------------------------------------------------------
try
    distX = give_Eulerian_Lagrangian_Distance(x(indX),xLag(ind_Lag),Lx);
    distY = give_Eulerian_Lagrangian_Distance(y(indY),yLag(ind_Lag'),Ly);
catch
    fprintf('\n\n\n - ERROR - \n');
    fprintf('\n\n - ERROR ERROR - \n');
    fprintf('\n\n - ERROR ERROR ERROR - \n');
    fprintf('\n\n - ERROR ERROR ERROR ERROR - \n\n\n');
    error('BLOW UP! (*forces TOO large*) -> try decreasing the time-step or decreasing material property stiffnesses');
end

%--------------------------------------------------------------------------------------
% Initialize delta_X, delta_Y matrices for storing delta-function info for each Lag Pt.
%--------------------------------------------------------------------------------------
delta_X = zeros(Nb, Nx);
delta_Y = zeros(Ny, Nb);

%-----------------------------------------------
% Get delta function kernels
%-----------------------------------------------
delta_1D_x = give_Delta_Kernel(distX, dx);
delta_1D_y = give_Delta_Kernel(distY, dy);


%-------------------------------------------------------
%               ORIGINAL IMPLEMENTATION
%  --> Tried case w/ 4096 Lag. Pts and it was faster
%      than more vectorized version below
%-------------------------------------------------------
[row,col] = size(ind_Lag);
for i=1:row
    for j=1:col
        
        % Get Eulerian/Lagrangian indices to use for saving non-zero delta-function values
        xID = indX(i,j);
        indy= ind_Lag(i,j);
        yID = indY(j,i);
        
        % Store non-zero delta-function values into delta_X / delta_Y matrices at correct indices!
        delta_X(indy,xID) = delta_1D_x(i,j);
        delta_Y(yID,indy) = delta_1D_y(j,i);
        
    end
end


%----------------------------------------------------
%       TRY TO GET AROUND NESTED FOR-LOOP ABOVE
%  --> Actually slower...
%  --> Tried case w/ 4096 Lag. Pts...
%----------------------------------------------------
% for i=1:Nb
%         
%         % Get Eulerian/Lagrangian indices to use for saving non-zero delta-function values
%         %xID = indX(i,:);
%         %indy= i;%ind_Lag(i,:)
%         %yID = indY(:,i);
%         
%         % Store non-zero delta-function values into delta_X / delta_Y matrices at correct indices!
%         delta_X(i, indX(i,:) ) = delta_1D_x(i,1:supp);
%         delta_Y( indY(:,i) ,i) = delta_1D_y(1:supp,i);
%         
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: CHECK if BEAM points have passed through boundary and then 
%           translates them appropriately for calculation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Xp_N,Xq_N,Xr_N] = check_If_Beam_Points_Pass_Through_Boundary(ds,Lx,Xp,Xq,Xr)

    % CHECKS FOR IF POINTS PASSED THRU BNDRY
    dX_pq = ( Xp - Xq );
    dX_qr = ( Xq - Xr );
    
    if abs(dX_pq) > 5*ds

        % MEANS point p has moved thru; check if moved through right/left bndry
        if dX_pq < 0
            Xp_N = Lx + Xp;
        else
            Xp_N = -Lx+Xp;
        end
        Xq_N = Xq;
        Xr_N = Xr;
       
    else
       Xp_N = Xp;
       Xq_N = Xq;
       Xr_N = Xr;
    end
    
    
    if abs(dX_qr) > 5*ds

        % MEANS point r has moved thru; check if moved through right/left bndry
        if dX_qr < 0
            Xr_N = -Lx+Xr;
        else
            Xr_N = Lx + Xr;
        end
        Xq_N = Xq;
        if abs(dX_pq) < 5*ds 
            Xp_N = Xp;
        end


    end

