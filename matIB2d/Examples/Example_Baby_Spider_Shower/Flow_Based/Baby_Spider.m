%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  battistn@tcnj.edu
% Date Created: May 27th, 2015
% Institution: UNC-CH
%
% Date Modified: 06/24/2021
% Institution Modified: TCNJ
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
% If you would like us to add a specific muscle model, please let Nick (battistn@tcnj.edu) know.
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates the BABY_SPIDER-EXAMPLE geometry and prints associated input files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Baby_Spider()

%
% Grid Parameters (MAKE SURE MATCHES IN input2d !!!)
%
Nx =  384;       % # of Eulerian Grid Pts. in x-Direction (MUST BE EVEN!!!)
Ny =  256;       % # of Eulerian Grid Pts. in y-Direction (MUST BE EVEN!!!)
Lx = 1.5;        % Length of Eulerian Grid in x-Direction
Ly = 1.0;        % Length of Eulerian Grid in y-Direction


% Immersed Structure Geometric / Dynamic Parameters %
Ls = 0.175;       % Length of baby-spider web
ds = Lx/(2*Nx);  % Lagrangian spacing, ds
struct_name = 'spider'; % Name for .vertex, .spring, etc files.


% Call function to construct geometry
[xLag,yLag,Ninfo] = give_Me_Immersed_Boundary_Geometry(ds,Ny,Lx,Ly,Ls);


% Plot Geometry to test
plot(xLag,yLag,'r-'); hold on;
plot(xLag,yLag,'*'); hold on;
xlabel('x'); ylabel('y');
axis([0 Lx 0 Ly]);


% Print INFO to screen
fprintf('\n\n                   INFO TO CHECK INPUT FILES!\n\n');
fprintf('# of Lag. Pts in WEB: %d\n',Ninfo(1));
fprintf('# of Lag. Pts in FLOOR: %d\n',Ninfo(2));
fprintf('Total # of Lag. Pts: %d\n',Ninfo(1)+Ninfo(2));
fprintf('Index of MASSIVE Pt (last pt. of WEB): %d\n\n\n',Ninfo(1));


% Prints .vertex file!
print_Lagrangian_Vertices(xLag,yLag,struct_name);


% Prints .mass file!
k_Mass = 1e4;        % 'spring' stiffness parameter for tethering
Mass =   0.0225;     % "MASS" value for 'MASSIVE' nodal movement
print_Lagrangian_Mass_Pts(xLag,k_Mass,Mass,struct_name,Ninfo);


% Prints .spring file!
k_Spring = 0.0375 / ds^2;  % k_{Non-Dim} =  0.037;
print_Lagrangian_Springs(xLag,yLag,k_Spring,struct_name,Ninfo,ds);


% Prints .beam file!
k_Beam = 0.38 / ds^2; % k_{Non-Dim} = 0.38;
C = 0.0;
print_Lagrangian_Beams(xLag,yLag,k_Beam,C,struct_name,Ninfo);


% Prints .target file! (for "Ground")
%k_Target = 2.5e5;
%print_Lagrangian_Target_Pts(xLag,k_Target,struct_name,Ninfo);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints VERTEX points to a file called spider.vertex
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Vertices(xLag,yLag,struct_name)

    N = length(xLag);

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
% FUNCTION: prints MASS points to a file called spider.mass
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
function print_Lagrangian_Mass_Pts(xLag,kMass,Mass,struct_name,Ninfo)

    %LOOP OVER LAG PTS FOR MASS PTS IN LAGRANGIAN INDEXING

    N = 1;  % ONLY 1 MASS PT.!

    mass_fid = fopen([struct_name '.mass'], 'w');

    fprintf(mass_fid, '%d\n', N );

    % Single Mass Pt. (bottom point of web!)
    fprintf(mass_fid, '%d %1.16f %1.16f\n', Ninfo(1), kMass,Mass);
    
    fclose(mass_fid); 
       
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints Vertex points to a file called spider.vertex
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Target_Pts(xLag,k_Target,struct_name,Ninfo)

    N = Ninfo(2);   % # of Lag. Pts on FLOOR!
    Ns= Ninfo(1)+1; % Starting Lag. Pt. INDEX for FLOOR
    
    target_fid = fopen([struct_name '.target'], 'w');

    fprintf(target_fid, '%d\n', N );

    %Loops over Lag. Pts on FLOOR ONLY
    for s = Ns:Ns+(N-1)
        fprintf(target_fid, '%d %1.16e\n', s, k_Target);
    end

    fclose(target_fid); 
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints BEAM (Torsional Spring) points to a file called spider.beam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Beams(xLag,yLag,k_Beam,C,struct_name,Ninfo)

    % k_Beam: beam stiffness
    % C: beam curvature
    
    N = Ninfo(1)-2; % NOTE: Total # of beams = % of Lag Pts. on WEB - 2

    beam_fid = fopen([struct_name '.beam'], 'w');

    fprintf(beam_fid, '%d\n', N );

    %BEAMS BETWEEN VERTICES
    for s = 2:N-1
        fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',s-1, s, s+1, k_Beam, C);  
    end
    fclose(beam_fid); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints SPRING points to a file called spider.spring
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Springs(xLag,yLag,k_Spring,struct_name,Ninfo,ds)

    N = Ninfo(1)-1;  % # of springs (just on 'web' of geometry)

    spring_fid = fopen([struct_name '.spring'], 'w');

    fprintf(spring_fid, '%d\n', N );

    %spring_force = kappa_spring*ds/(ds^2);

    %SPRINGS BETWEEN VERTICES ON WEBBING ONLY
    for s = 1:N
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, s+1, k_Spring, ds);  
    end
    fclose(spring_fid); 
    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates the Lagrangian structure geometry
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag,Ninfo] = give_Me_Immersed_Boundary_Geometry(ds,Ny,Lx,Ly,Ls)


%
% Create Spider Web Geometry (springs, beams, and 1-mass pt)
%
yS =   0.85*Ly;                       % Highest Point of Web
yWeb = yS:-ds:yS-Ls;                  % yPts of web
Nweb = length(yWeb);                  % # of Lag. Pts in Web
xWeb = 0.25*ones(1,Nweb);           % xPts of web

%
% Create Floor (for target pts)
% 
frac = 0.025;                                     % Fraction of bottom of domain, 'not floored'
frac2 = 0.01;                                     % Fraction of height of domain for floor to be placed
xFloor = (frac/2)*Lx:2*ds:Lx - (frac/2);  % xPts of floor
Nfloor = length(xFloor);                          % # of Lag. Pts in Floor
yFloor = frac2*Ly*ones(1,Nfloor);                 % yPts of floor

%
% NO "GROUND" (NO TARGET POINTS)
%
xLag = [xWeb];             %Vector of all x-Lag. Pts.
yLag = [yWeb];             %Vector of all y-Lag. Pts.
Ninfo= [Nweb 0 Nweb]; %Vector of # of Lag. Pts for each part of geometry

%
% INCLUDE "GROUND" (TARGET POINTS)
%
%xLag = [xWeb xFloor];             %Vector of all x-Lag. Pts.
%yLag = [yWeb yFloor];             %Vector of all y-Lag. Pts.
%Ninfo= [Nweb Nfloor Nweb+Nfloor]; %Vector of # of Lag. Pts for each part of geometry

