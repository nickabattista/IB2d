%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  nick.battista@unc.edu
% Date Created: September 9th, 2016
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates the geometry and prints associated input files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Make_Your_Geometry_and_Input_Files()

%
% Grid Parameters (MAKE SURE MATCHES IN input2d !!!)
%
Nx =  1024;      % # of Eulerian Grid Pts. in x-Direction (MUST BE EVEN!!!)
Lx = 5.0;        % Length of Eulerian Grid in x-Direction
dx = Lx/Nx;      % Grid spatial resolution

%
% Immersed Structure Geometric / Dynamic Parameters %
%
ds = 0.45*dx;               % Lagrangian Pt. Spacing (2x resolution of Eulerian grid)
struct_name = 'insect_HT'; % Name for .vertex, .spring, etc files. (must match what's in 'input2d')


% Call function to construct geometry
[xLag,yLag] = give_Me_Immsersed_Boundary_Geometry(ds,Nx,Lx);
yLag = yLag - 2+0.125;

% Plot Geometry to test
plot(xLag,yLag,'r-'); hold on;
plot(xLag,yLag,'*'); hold on;
xlabel('x'); ylabel('y');
axis([0 Lx 0 Lx]);



% Prints .vertex file!
print_Lagrangian_Vertices(xLag,yLag,struct_name);


% Prints .spring file!
%k_Spring = 2.5e4;                    % Spring stiffness (does not need to be equal for all springs)
%ds_Rest = 0.0;                       % Spring resting length (does not need to be equal for all springs)
%print_Lagrangian_Springs(xLag,yLag,k_Spring,ds_Rest,struct_name);


% Prints .beam file!
% k_Beam = 0.5;                      % Beam Stiffness (does not need to be equal for all beams)
% C = compute_Curvatures(xLag,yLag)  % Computes curvature of initial configuration
%print_Lagrangian_Beams(xLag,yLag,k_Beam,C,struct_name);


% Prints .target file!
k_Target = 2e5;
print_Lagrangian_Target_Pts(xLag,k_Target,struct_name);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints VERTEX points to a file called 'struct_name'.vertex
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Vertices(xLag,yLag,struct_name)

    N = length(xLag); % Total # of Lag. Pts

    vertex_fid = fopen([struct_name '.vertex'], 'w');

    fprintf(vertex_fid, '%d\n', N );

    %Loops over all Lagrangian Pts.
    for s = 1:N
        X_v = xLag(s);
        Y_v = yLag(s);
        fprintf(vertex_fid, '%1.16e %1.16e\n', X_v, Y_v);
    end

    fclose(vertex_fid);

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints TARGET points to a file called 'struct_name'.target
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Target_Pts(xLag,k_Target,struct_name)

    N = length(xLag);

    target_fid = fopen([struct_name '.target'], 'w');

    fprintf(target_fid, '%d\n', N );

    %Loops over all Lagrangian Pts.
    for s = 1:N
        fprintf(target_fid, '%d %1.16e\n', s, k_Target);
    end

    fclose(target_fid);   
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints SPRING points to a file called 'struct_name'.spring
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Springs(xLag,yLag,k_Spring,ds_Rest,struct_name)

    N = length(xLag);

    spring_fid = fopen([struct_name '.spring'], 'w');

    fprintf(spring_fid, '%d\n', N );    % Print # of springs 

    %SPRINGS BETWEEN VERTICES
    for s = 1:N
            if s < N         
                fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, s+1, k_Spring, ds_Rest);  
            else
                %Case s=N
                fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, 1,   k_Spring, ds_Rest);  
            end
    end
    fclose(spring_fid);      

    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints BEAM (Torsional Spring) points to a file called 'struct_name'.beam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Beams(xLag,yLag,k_Beam,C,struct_name)

    % k_Beam: beam stiffness
    % C: beam curvature
    
    N = length(xLag); % NOTE: Total number of beams = Number of Total Lag Pts. - 2

    beam_fid = fopen([struct_name '.beam'], 'w');

    fprintf(beam_fid, '%d\n', N );

    %spring_force = kappa_spring*ds/(ds^2);

    %BEAMS BETWEEN VERTICES
    for s = 2:N-1
            if  s <= N-1         
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',s-1, s, s+1, k_Beam, C(s) );  
            else
                %Case s=N
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',s-1, s, 1,   k_Beam, C(s) );  
            end
    end
    fclose(beam_fid); 
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes "curvature" of starting configuration
% 
% NOTE: not curvature in the traditional geometric sense, in the 'discrete'
% sense through cross product.
%
% NOTE: assumes a CLOSED structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C = compute_Curvatures(xLag,yLag)

N = length(xLag);
C = zeros( N );

%Note: needs to be done same order as you print .beam file!
for i=1:N
   
   % Pts Xp -> Xq -> Xr (same as beam force calc.)
   
   if ( (i > 1) && (i < N) )
   
        Xp = xLag(i-1); Xq = xLag(i); Xr = xLag(i+1);
        Yp = yLag(i-1); Yq = yLag(i); Yr = yLag(i+1);
   
   elseif (i==1)
       
        Xp = xLag(N); Xq = xLag(i); Xr = xLag(i+1);
        Yp = yLag(N); Yq = yLag(i); Yr = yLag(i+1);
       
   elseif (i==N)
       
        Xp = xLag(N-1); Xq = xLag(N); Xr = xLag(1);
        Yp = yLag(N-1); Yq = yLag(N); Yr = yLag(1);
       
   end
       
   C(i) = (Xr-Xq)*(Yq-Yp) - (Yr-Yq)*(Xq-Xp); %Cross product btwn vectors
      
end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates the Lagrangian structure geometry
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [xLag,yLag] = give_Me_Immsersed_Boundary_Geometry(ds,Nx,Lx)
 
% ds: Lagrangian pt. spacing
% Nx: Eulerian grid resolution

% Gives elliptical subset for tube side
a = 0.5;    % semi-major axis (horizontal)
b = 0.45;   % semi-minor axis (vertical) 
ang = pi/6; % angle start/cuttoff for ellipse
[xE,yE] = give_Me_Elliptical_Geometry(ds,a,b,ang);
xOff = 2*xE(1);

% Gives top and bottom portion of one section of tube
d = 0.2;                  % separation distance
xTube = [xE+3*xOff xE+3*xOff xE+2*xOff xE+2*xOff xE+xOff xE+xOff xE xE];          % top first, then bottom
xTube = xTube + 1;
xTube = xTube(end:-1:1);  % reverse it!
yTube = [yE+d/2 -yE-d/2 yE+d/2 -yE-d/2 yE+d/2 -yE-d/2 yE+d/2 -yE-d/2]; % top first, then bottom 
yTube = yTube + Lx/2;

% Give indices for to attach valves
ind_Top = 1;                 % node index for first point on top (on RHS)
ind_Bot = length(xTube)/8+1; % node index for first point on bottom (on RHS)

xFlat = xTube(ind_Top)-Lx/20:ds:xTube(ind_Top); 
yFlatTop = yTube(ind_Top)*ones(1,length(xFlat));
yFlatBot = yTube(ind_Bot)*ones(1,length(xFlat)); 

xLag = [xTube xFlat xFlat];
yLag = [yTube yFlatTop yFlatBot];

% Gives valve geometry
[xV,yV] = give_Valve_Geometry(d,ds);

% TESTING (i.e., PLOT IT, yo!)
%plot(xTube,yTube,'b*'); hold on
%len = length(xTube);
%plot(xTube(len/8+1),yTube(len/8+1),'r*'); hold on
%plot(xTube(1),yTube(1),'k*'); hold on;
%axis([0 5 0 5]);
%pause();

%plot(xV,yV,'g*'); hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives valve geometry inside a distance, d, between the edges of
% the tube
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xV,yV] = give_Valve_Geometry(d,ds)

% for A*x^2 geometry -> 2 eqns and 2 unknowns
%y(0) = 0 = 0
%y(0.9d/2) = A(0.9/2)^2*d^2 = 1.5*d
%
% ----> A = 1.5 / ( (0.45)^2 * d )
%
A = 1.0 / ( (0.45)^2 * d );

yV = -0:ds/2:0.9/2*d;
xV = -A*yV.^2;


%figure(2)
%plot(xV,yV,'b*'); hold on;
%pause();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives elliptical geometry
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xE,yE] = give_Me_Elliptical_Geometry(ds,a,b,ang)

dtheta = ds/max(a,b); % use minimum of axis lengths

n = 1;
xE(n) = a*cos(ang); yE(n) = b*sin(ang);
ang_n = ang + dtheta;

while ang_n < (pi/2)
    
    n = n+1;                % update counter
    xE(n) = a*cos(ang_n);   % next xVal along ellipse
    yE(n) = b*sin(ang_n);   % next yVal along ellipse
    ang_n = ang_n + dtheta; % update angle
    
end
xE(n+1) = a*cos(pi/2);
yE(n+1) = b*sin(pi/2);

xE = [xE -xE(end-1:-1:1)];
yE = [yE yE(end-1:-1:1)];

yE = yE - yE(1);           % Translate down

%plot(xE,yE,'*'); hold on;
%axis([-a a -a a]);

xLength = abs( 2*xE(1) );
fprintf('\nThe x-Length of the tube section is: %d\n',xLength);


