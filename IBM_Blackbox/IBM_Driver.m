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
pDump= grid_Info(8); % Print (Plot) Dump interval


% MODEL STRUCTURE DATA STORED %
springs_Yes = model_Info(1);         % Springs: 0 (for no) or 1 (for yes) 
update_Springs_Flag = model_Info(2); % Update_Springs: 0 (for no) or 1 (for yes)
target_pts_Yes = model_Info(3);      % Target_Pts: 0 (for no) or 1 (for yes)
update_Target_Pts = model_Info(4);   % Update_Target_Pts: 0 (for no) or 1 (for yes)
beams_Yes = model_Info(5);           % Beams: 0 (for no) or 1 (for yes)
update_Beams_Flag = model_Info(6);   % Update_Beams: 0 (for no) or 1 (for yes)
muscles_Yes = model_Info(7);         % Muscles: 0 (for no) or 1 (for yes)
arb_ext_force_Yes = model_Info(8);   % Arbitrary External Force: 0 (for no) or 1 (for yes)

%Lagrangian Structure Data
ds = Lx / (2*Nx);                   %Lagrangian Spacing
grid_Info(9) = ds;


% Create EULERIAN Mesh (these assume periodicity in x and y)
x = (0:dx:Lx-dx);  X = [];
y = (0:dy:Ly-dy)'; Y = [];
%Create x-Mesh
for i=1:Nx
    X = [X; x]; 
end
%Create y-Mesh
for i=1:Ny
    Y = [Y y];
end




% % % % % HOPEFULLY WHERE I CAN READ IN INFO!!! % % % % %


% READ IN LAGRANGIAN POINTS %
[Nb,xLag,yLag] = read_Vertex_Points(struct_name);
grid_Info(8) = Nb;          % # Total Number of Lagrangian Pts.
xLag_P = xLag;              % Initialize previous Lagrangian x-Values (for use in muscle-model)
yLag_P = yLag;              % Initialize previous Lagrangian y-Values (for use in muscle-model)


% READ IN SPRINGS (IF THERE ARE SPRINGS) %
if ( springs_Yes == 1 )
    springs_info = read_Spring_Points(struct_name);
        %springs_info: col 1: starting spring pt (by lag. discretization)
        %              col 2: ending spring pt. (by lag. discretization)
        %              col 3: spring stiffness
        %              col 4: spring resting lengths
else
    springs_info = 0;  %just to pass placeholder into "please_Find_Lagrangian_Forces_On_Eulerian_grid function"
end




% READ IN MUSCLES (IF THERE ARE MUSCLES) %
if ( muscles_Yes == 1 )
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





% READ IN TARGET POINTS (IF THERE ARE TARGET PTS) %
if ( target_pts_Yes == 1)
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




% READ IN BEAMS (IF THERE ARE BEAMS) %
if ( beams_Yes == 1)
    beams_info = read_Beam_Points(struct_name);
    %beams:      col 1: 1ST PT.
    %            col 2: MIDDLE PT. (where force is exerted)
    %            col 3: 3RD PT.
    %            col 4: beam stiffness
    %            col 5: curavture
else
    beams_info = 0;
end




% Initialize the initial velocities to zero.
U = zeros(Ny,Nx);
V = U;


% ACTUAL TIME-STEPPING IBM SCHEME!
ct = 0;
while current_time < T_FINAL
    
    % Step 1: Update Position of Boundary of membrane at half time-step
    % Variables end with h if it is a half-step
    [xLag_h, yLag_h] = please_Move_Lagrangian_Point_Positions(U, V, xLag, yLag, xLag, yLag, x, y, dt/2, grid_Info);
    
    if ( ( update_Springs_Flag == 1 ) && ( springs_Yes == 1 ) )
       springs_info = update_Springs(dt,current_time,xLag,springs_info); 
    end
    
    if ( ( update_Target_Pts == 1 ) && ( target_pts_Yes == 1) )
       target_info = update_Target_Point_Positions(dt,current_time,target_info); 
    end
    
    if ( ( update_Beams_Flag == 1 ) && ( Beams_Yes == 1) )
       beams_info = update_Beams(dt,current_time,beams_info); 
    end
    
    
    % Step 2: Calculate Force coming from membrane at half time-step
    [Fxh, Fyh] =           please_Find_Lagrangian_Forces_On_Eulerian_grid(dt, current_time, xLag_h, yLag_h, xLag_P, yLag_P, x, y, grid_Info, model_Info, springs_info, target_info, beams_info, muscles_info);
    
    if arb_ext_force_Yes == 1 
        [Fx_Arb, Fy_Arb] = please_Compute_External_Forcing(dt, current_time, x, y, grid_Info, U, V);
        Fxh = Fxh + Fx_Arb;
        Fyh = Fyh + Fy_Arb;
    end
    
    
    % Step 3: Solve for Fluid motion
    [Uh, Vh, U, V] =   please_Update_Fluid_Velocity(U, V, Fxh, Fyh, rho, mu, grid_Info, dt);

    % Step 4: Update Position of Boundary of membrane again for a half time-step
    xLag_P = xLag_h;     % Stores old Lagrangian x-Values (for muscle model)
    yLag_P = yLag_h;     % Stores old Lagrangian y-Values (for muscle model)
    [xLag, yLag] =     please_Move_Lagrangian_Point_Positions(Uh, Vh, xLag, yLag, xLag_h, yLag_h, x, y, dt, grid_Info);

    % Plot Lagrangian/Eulerian Dynamics!
    if mod(ct,pDump) ==0
        vort = give_Me_Vorticity(U,V,dx,dy);
        uMag = give_Me_Magnitude_Velocity(U,V);
        please_Plot_Results(dx,dy,X,Y,U,V,vort,xLag,yLag);
    end

    
    % Update current_time value
    current_time = current_time+dt;
    fprintf('Current Time(s): %6.6f\n',current_time);
    
    %ct = ct + 1;
    %if mod(ct,200)==0
    %    pause();
    %end
    
end %ENDS TIME-STEPPING LOOP






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
% FUNCTION: Reads in the # of springs and all MASTER NODEs, SLAVE NODEs,
%           spring STIFFNESSES, spring RESTING LENGTHS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function springs = read_Spring_Points(struct_name)

filename = [struct_name '.spring'];  %Name of file to read in

fileID = fopen(filename);

    % Read in the file, use 'CollectOutput' to gather all similar data together
    % and 'CommentStyle' to to end and be able to skip lines in file.
    C = textscan(fileID,'%f %f %f %f','CollectOutput',1);

fclose(fileID);        %Close the data file.

spring_info = C{1};    %Stores all read in data in vertices (N+1,2) array

%Store all elements on .spring file into a matrix starting w/ 2nd row of read in data.
springs = spring_info(2:end,1:4);

%springs: col 1: starting spring pt (by lag. discretization)
%         col 2: ending spring pt. (by lag. discretization)
%         col 3: spring stiffness
%         col 4: spring resting lengths




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
    
    
    