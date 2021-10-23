%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  battistn@tcnj.edu
% Date Simulation Created: August 28, 2019
% Institution: TCNJ
% 
% Date IB2d Was Created: May 27th, 2015
% Institution: UNC-CH
%
% This code is capable of creating Lagrangian Structures using:
% 	1. Springs
% 	2. Beams (*torsional springs)
% 	3. Target Points
%	4. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-Tension)")
%   .
%   .
%   .
%
% One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting lengths, etc
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
% 
% If you would like us to add a specific fiber model or example, please let Nick (battistn@tcnj.edu) know.
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates the PENDULUM POINT MASS EXAMPLE geometry and prints associated input files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Give_Pendulum_Geometry_PointMass(Mass)

%
% Grid Parameters (MAKE SURE MATCHES IN input2d !!!)
%
Nx =  1024;      % # of Eulerian Grid Pts. in x-Direction (MUST BE EVEN!!!)
Ny =  1024;      % # of Eulerian Grid Pts. in y-Direction (MUST BE EVEN!!!)
Lx = 1.0;        % Length of Eulerian Grid in x-Direction
Ly = 1.0;        % Length of Eulerian Grid in y-Direction


% Immersed Structure Geometric / Dynamic Parameters %
ds= min(Lx/(2*Nx),Ly/(2*Ny));  % Lagrangian spacing
Lp = 0.2;                      % Length of Pendulum                
struct_name = 'pendulum';      % Name for .vertex, .spring, etc files.

% x and y vectors represent pairs of coordinates describing pendulum and its point mass bob.
x = -Lp/2:Lp:Lp/2;
y = zeros(length(x));

% Angle of rotation from x-axis
t = pi/5;

% Carries out rotation
for j=1:length(x)
    xN(j) = x(j)*cos(t) - y(j)*sin(t);
    yN(j) = x(j)*sin(t) + y(j)*cos(t);    
end

%
% Shifts to appropriate place in domain
%
xN = xN - xN(end) + Lx/2;
yN = yN - yN(end) + 2*Ly/3;

%
% Plots the Geometry
%
plot(x, y, 'b.-'); hold on;
plot(xN, yN, 'r.-'); hold on;
axis([0 1 0 1]);

%
% Defines xLag and yLag for printing input files below
%
xLag = xN;
yLag = yN;

% Prints .vertex file!
print_Lagrangian_Vertices(xLag,yLag,struct_name);

% Prints .mass file!
k_Mass = 2.5e6;            % 'spring' stiffness parameter for mass tethering
%Mass = 1e2;               % "MASS" value for point mass (input to this script)
print_Lagrangian_Mass_Pts(xLag,k_Mass,Mass,struct_name);

% Prints .spring file!
k_Spring = 5e7; ds_Rest = sqrt( ( xLag(1)-xLag(2)  )^2 + ( yLag(1)-yLag(2)  )^2    );
print_Lagrangian_Springs(xLag,yLag,k_Spring,ds_Rest,struct_name);

% Prints .target file!
k_Target = 5e7;
print_Lagrangian_Target_Pts(xLag,k_Target,struct_name);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints VERTEX points to a file called rubberband.vertex
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
% FUNCTION: prints TARGET points to a file called rubberband.vertex
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Target_Pts(xLag,k_Target,struct_name)

    N = 1; %for 3 channels

    target_fid = fopen([struct_name '.target'], 'w');

    fprintf(target_fid, '%d\n', N );

    %Loops over all Lagrangian Pts.
   
    fprintf(target_fid, '%d %1.16e\n', 2, k_Target); %to choose correct lag. pt.
 

    fclose(target_fid); 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints TARGET points to a file called struct_name.mass
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
function print_Lagrangian_Mass_Pts(xLag,kMass,Mass,struct_name)

    %LOOP OVER LAG PTS FOR MASS PTS IN LAGRANGIAN INDEXING

    N = 1;

    mass_fid = fopen([struct_name '.mass'], 'w');

    fprintf(mass_fid, '%d\n', N );

    %Loops over all Lagrangian Pts.
   
    fprintf(mass_fid, '%d %1.16e %1.16e\n', 1, kMass, Mass);
    
    fclose(mass_fid); 
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints BEAM (Torsional Spring) points to a file called rubberband.beam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Beams(xLag,yLag,k_Beam,C,struct_name)

    % k_Beam: beam stiffness
    % C: beam curvature
    
    N = length(xLag); % NOTE: Total number of beams = Number of Total Lag Pts. - 2

    beam_fid = fopen([struct_name '.beam'], 'w');

    fprintf(beam_fid, '%d\n', 3*N );

    %spring_force = kappa_spring*ds/(ds^2);

    %BEAMS BETWEEN VERTICES
    for s = 1:N
            if s==1 
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',N, s, s+1,   k_Beam, C);  
            elseif  s <= N-1         
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',s-1, s, s+1, k_Beam, C);  
            elseif s==N
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',s-1, s, 1,   k_Beam, C);
            end
    end
    
    for s = N+1:2*N
            if s==N+1 
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',2*N, s, s+1,   k_Beam, C);  
            elseif  s <= 2*N-1         
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',s-1, s, s+1, k_Beam, C);  
            elseif s==2*N
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',s-1, s, N+1,   k_Beam, C);
            end
    end
    
    for s = 2*N+1:3*N
            if s==2*N+1 
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',3*N, s, s+1,   k_Beam, C);  
            elseif  s <= 3*N-1         
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',s-1, s, s+1, k_Beam, C);  
            elseif s==3*N
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',s-1, s, 2*N+1,   k_Beam, C);
            end
    end
    fclose(beam_fid); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints SPRING points to a file called rubberband.spring
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Springs(xLag,yLag,k_Spring,ds_Rest,struct_name)

    N = length(xLag);

    spring_fid = fopen([struct_name '.spring'], 'w');

    fprintf(spring_fid, '%d\n', 1 ); %N MUST BE EVEN!

    %spring_force = kappa_spring*ds/(ds^2);

    %SPRINGS BETWEEN VERTICES
                 
    fprintf(spring_fid, '%d %d %1.16e %1.16e\n', 1, 2, 2*k_Spring, ds_Rest);  

    fclose(spring_fid); 
    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates the Lagrangian structure geometry for cylinder
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag] = give_Me_Circular_Immersed_Boundary_Geometry(ds,r,x0,y0)

% The immsersed structure is a circle %

dtheta = ds/ (2*r);
theta = 0; i=1;
while theta < 2*pi
   xLag(i) = x0 - r*cos(theta);
   yLag(i) = y0 - r*sin(theta);
   theta = theta + dtheta;
   i=i+1;
end

   
   
