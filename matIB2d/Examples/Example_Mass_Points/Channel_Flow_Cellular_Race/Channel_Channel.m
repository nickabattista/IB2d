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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates the CHANNEL_CHANNEL-EXAMPLE geometry and prints associated input files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Channel_Channel()

%
% Grid Parameters (MAKE SURE MATCHES IN input2d !!!)
%
Nx =  64;        % # of Eulerian Grid Pts. in x-Direction (MUST BE EVEN!!!)
Ny =  64;        % # of Eulerian Grid Pts. in y-Direction (MUST BE EVEN!!!)
Lx = 1.0;        % Length of Eulerian Grid in x-Direction
Ly = 1.0;        % Length of Eulerian Grid in y-Direction


% Immersed Structure Geometric / Dynamic Parameters %
ds= min(Lx/(2*Nx),Ly/(2*Ny));  % Lagrangian spacing
L = 0.9*Lx;                    % Length of Channel
w = 0.2*Ly;                    % Width of Channel
x0 = 0.3;                      % x-Center for Cylinder
y0 = 0.5;                      % y-Center for Cylinder
r = w/6;                       % Radii of Cylinder
struct_name = 'channel'; % Name for .vertex, .spring, etc files.


% Call function to construct geometry
[xLag_T,yLag_T] = give_Me_Channel_Immsersed_Boundary_Geometry(ds,L,w,Lx,Ly);
[xLag_C,yLag_C] = give_Me_Cylinder_Immsersed_Boundary_Geometry(ds,r,x0,y0);

%Upper Channel / Cell
xLag_T1 = xLag_T; yLag_T1 = yLag_T + 0.3;
xLag_C1 = xLag_C; yLag_C1 = yLag_C + 0.3;

%Middle Channel / Cell
xLag_T2 = xLag_T; yLag_T2 = yLag_T;
xLag_C2 = xLag_C; yLag_C2 = yLag_C;

%Lower Channel / Cell
xLag_T3 = xLag_T; yLag_T3 = yLag_T - 0.3;
xLag_C3 = xLag_C; yLag_C3 = yLag_C - 0.3;


% Plot Geometry to test
plot(xLag_T1(1:end/2),yLag_T1(1:end/2),'r-'); hold on;
plot(xLag_T2(1:end/2),yLag_T2(1:end/2),'r-'); hold on;
plot(xLag_T3(1:end/2),yLag_T3(1:end/2),'r-'); hold on;
plot(xLag_T1(end/2+1:end),yLag_T1(end/2+1:end),'r-'); hold on;
plot(xLag_T2(end/2+1:end),yLag_T2(end/2+1:end),'r-'); hold on;
plot(xLag_T3(end/2+1:end),yLag_T3(end/2+1:end),'r-'); hold on;
%
plot(xLag_C1,yLag_C1,'r-'); hold on;
plot(xLag_C2,yLag_C2,'r-'); hold on;
plot(xLag_C3,yLag_C3,'r-'); hold on;
%
plot(xLag_T1,yLag_T1,'*'); hold on;
plot(xLag_T2,yLag_T2,'*'); hold on;
plot(xLag_T3,yLag_T3,'*'); hold on;
%
plot(xLag_C1,yLag_C1,'g*'); hold on;
plot(xLag_C2,yLag_C2,'g*'); hold on;
plot(xLag_C3,yLag_C3,'g*'); hold on;
%
xlabel('x'); ylabel('y');
axis([0 Lx 0 Ly]);


% Combine all lagrangian pts into one vector
xLag = [xLag_C1 xLag_C2 xLag_C3 xLag_T1 xLag_T2 xLag_T3]; % Add xLagPts from Circle to xLag Pt. Vector (*no springs or beams*)
yLag = [yLag_C1 yLag_C2 yLag_C3 yLag_T1 yLag_T2 yLag_T3]; % Add xLagPts from Circle to xLag Pt. Vector (*no springs or beams*)

% Prints .vertex file!
print_Lagrangian_Vertices(xLag,yLag,struct_name);

% Prints .mass file!
k_Mass = 1e6;         % 'spring' stiffness parameter for tethering
Mass = [1e-2 5e-1 1e0];          % "MASS" value for 'ghost' nodal movement
print_Lagrangian_Mass_Pts(xLag_C1,k_Mass,Mass,struct_name);

% Prints .spring file!
k_Spring = 5e7; ds_Rest = sqrt( ( xLag_C(1)-xLag_C(2)  )^2 + ( yLag_C(1)-yLag_C(2)  )^2    );
print_Lagrangian_Springs(xLag_C,yLag_C,k_Spring,ds_Rest,r,struct_name);

% Prints .beam file!
k_Beam = 1e7; C = 0;%1/r; %(curvature of a circle)
print_Lagrangian_Beams(xLag_C,yLag_C,k_Beam,C,struct_name);

% Prints .target file!
k_Target = 1e7; Noff = 3*length(xLag_C);
print_Lagrangian_Target_Pts(xLag_T,Noff,k_Target,struct_name);

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

function print_Lagrangian_Target_Pts(xLag,Noff,k_Target,struct_name)

    N = 3*length(xLag); %for 3 channels

    target_fid = fopen([struct_name '.target'], 'w');

    fprintf(target_fid, '%d\n', N );

    %Loops over all Lagrangian Pts.
    for s = 1:N
        fprintf(target_fid, '%d %1.16e\n', Noff+s, k_Target); %to choose correct lag. pt.
    end

    fclose(target_fid); 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints TARGET points to a file called struct_name.mass
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
function print_Lagrangian_Mass_Pts(xLag,kMass,Mass,struct_name)

    %LOOP OVER LAG PTS FOR MASS PTS IN LAGRANGIAN INDEXING

    N = length(xLag);

    mass_fid = fopen([struct_name '.mass'], 'w');

    fprintf(mass_fid, '%d\n', 3*N );

    %Loops over all Lagrangian Pts.
    for s = 1:N
        fprintf(mass_fid, '%d %1.16e %1.16e\n', s, kMass,Mass(1));
    end
    for s = N+1:2*N
        fprintf(mass_fid, '%d %1.16e %1.16e\n', s, kMass,Mass(2));
    end
    for s = 2*N+1:3*N
        fprintf(mass_fid, '%d %1.16e %1.16e\n', s, kMass,Mass(3));
    end
    

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

function print_Lagrangian_Springs(xLag,yLag,k_Spring,ds_Rest,r,struct_name)

    N = length(xLag);

    spring_fid = fopen([struct_name '.spring'], 'w');

    fprintf(spring_fid, '%d\n', 3*N + 3/2*N ); %N MUST BE EVEN!

    %spring_force = kappa_spring*ds/(ds^2);

    %SPRINGS BETWEEN VERTICES
    for s = 1:3*N
            if s < N         
                fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, s+1, k_Spring, ds_Rest);  
            elseif s==N
                %Case s=N
                fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, 1,   k_Spring, ds_Rest);  
            elseif ( ( s>N ) && ( s<2*N ) )
                fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, s+1, k_Spring, ds_Rest);  
            elseif s==2*N
                fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, N+1,   k_Spring, ds_Rest); 
            elseif ( ( s>2*N ) && ( s<3*N ) )
                fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, s+1, k_Spring, ds_Rest);  
            elseif s==3*N
                fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, 2*N+1,   k_Spring, ds_Rest); 
            end
    end
    
    for s=1:N/2
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, s+N/2, k_Spring/100, 2*r);  
    end
    for s=1:N/2
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', N+s, N + s+N/2, k_Spring/100, 2*r);  
    end
    for s=1:N/2
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', 2*N+s, 2*N + s+N/2, k_Spring/100, 2*r);  
    end
    fclose(spring_fid); 
    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates the Lagrangian structure geometry for channel
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag] = give_Me_Channel_Immsersed_Boundary_Geometry(ds,L,w,Lx,Ly)

% The immsersed structure is a channel %
x = (Lx-L)/2:ds:(L+(Lx-L)/2);  %xPts
yBot = (Ly-w)/2;               %yVal for bottom of Channel
yTop = Ly - (Ly-w)/2;          %yVal for top of Channel

xLag = [x x];
yLag = [yBot*ones(1,length(x)) yTop*ones(1,length(x))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates the Lagrangian structure geometry for cylinder
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag] = give_Me_Cylinder_Immsersed_Boundary_Geometry(ds,r,x0,y0)

% The immsersed structure is a cylinder %

dtheta = ds/ (2*r);
theta = 0; i=1;
while theta < 2*pi
   xLag(i) = x0 - r*cos(theta);
   yLag(i) = y0 - r*sin(theta);
   theta = theta + dtheta;
   i=i+1;
end

   
   
