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

function Make_Your_Plates_and_Input_Files()

%
% Grid Parameters (MAKE SURE MATCHES IN input2d !!!)
%
Nx =  64;       % # of Eulerian Grid Pts. in x-Direction (MUST BE EVEN!!!)
Lx = 1.0;        % Length of Eulerian Grid in x-Direction
dx = Lx/Nx;      % Grid spatial resolution

%
% Immersed Structure Geometric / Dynamic Parameters %
%
ds = 0.5*dx;              % Lagrangian Pt. Spacing (2x resolution of Eulerian grid)
struct_name = 'corals';   % Name for .vertex, .spring, etc files. (must match what's in 'input2d')


% Call function to construct geometry
[xLag,yLag] = give_Me_Immsersed_Boundary_Geometry(ds,Nx,Lx);


% Plot Geometry to test
plot(xLag,yLag,'r-'); hold on;
plot(xLag,yLag,'*'); hold on;
xlabel('x'); ylabel('y');
axis([0 Lx 0 Lx]);


% Prints .vertex file!
print_Lagrangian_Vertices(xLag,yLag,struct_name);


% Prints .spring file!
% k_Spring = 2.5e4;                    % Spring stiffness (does not need to be equal for all springs)
% ds_Plates = dist;                    % Spring resting length (does not need to be equal for all springs)
% print_Lagrangian_Springs(xLag,yLag,k_Spring,ds_Plates,ds,struct_name);


% Prints .beam file!
% k_Beam = 0.5;                      % Beam Stiffness (does not need to be equal for all beams)
% C = compute_Curvatures(xLag,yLag)  % Computes curvature of initial configuration
%print_Lagrangian_Beams(xLag,yLag,k_Beam,C,struct_name);


% Prints .target file!
k_Target = 1e7;
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
% FUNCTION: prints SPRING points to a file called 'struct_name'.spring
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Springs(xLag,yLag,k_Spring,ds_Plates,ds,struct_name)

    N = length(xLag);

    spring_fid = fopen([struct_name '.spring'], 'w');

    fprintf(spring_fid, '%d\n', N-2 + N/2 );    % Print # of springs 

    %SPRINGS BETWEEN VERTICES
    for s = 1:N
            if s < N/2         
                fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, s+1, k_Spring, ds);  
            elseif s > N/2
                fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, s+1, k_Spring, ds);  
            end
    end
    
    % SPRINGS ACROSS PLATES
    for s=1:N/2
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, s+N/2, k_Spring, ds_Plates);  
    end
    
    fclose(spring_fid);     
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints TARGET points to a file called 'struct_name'.vertex
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
 
function [xLag,yLag] = give_Me_Immsersed_Boundary_Geometry(ds,Nx,L)

% ds: Lagrangian pt. spacing
% Nx: Eulerian grid resolution
% L:  Length of computational domain

ratio = 1024 / (2*Nx);

% Gives LEFT Coral
[x1,y1,N] = please_Give_Single_Polyp_Geometry(ratio);
x1 = x1 + L/10;

% Makes RIGHT Coral
x2 = x1 + L/5;
y2 = y1;

% TESTS GEOMETRY BY PLOTTING
plot(x1,y1,'k*'); hold on;
plot(x1(1),y1(1),'g*'); hold on;
plot(x1(end/2),y1(end/2),'b*'); hold on;
plot(x1(end/2+1),y1(end/2+1),'g*'); hold on;
plot(x1(end),y1(end),'b*'); hold on;
%
plot(x2,y2,'r*'); hold on;
plot(x2(1),y2(1),'g*'); hold on;
plot(x2(end/2),y2(end/2),'b*'); hold on;
plot(x2(end/2+1),y2(end/2+1),'g*'); hold on;
plot(x2(end),y2(end),'b*'); hold on;

xRef = [x1 x2]; yRef = [y1 y2]; 
ang = pi/4;

% Store Values for Centers of Rotation
xL_1 = xRef(1);       yL_1 = yRef(1);           % Left side of Left Pair
xR_1 = xRef(end/4+1); yR_1 = yRef(end/4+1);     % Right side of Left Pair

xL_2 = xRef(end/2+1);   yL_2 = yRef(end/2+1);   % Left side of Right Pair
xR_2 = xRef(3*end/4+1); yR_2 = yRef(3*end/4+1); % Right side of Right Pair

% -> Rotate Geometry <- %
% LEFT PAIR %
[xR_Ref_1,yR_Ref_1] = rotate_Geometry(ang,xR_1,yR_1,xRef(end/4+1:end/2),yRef(end/4+1:end/2) );
[xL_Ref_1,yL_Ref_1] = rotate_Geometry(-ang,xL_1,yL_1,xRef(1:end/4),yRef(1:end/4) );
% RIGHT PAIR %
[xR_Ref_2,yR_Ref_2] = rotate_Geometry(ang,xR_2,yR_2,xRef(3*end/4+1:end),yRef(3*end/4+1:end) );
[xL_Ref_2,yL_Ref_2] = rotate_Geometry(-ang,xL_2,yL_2,xRef(end/2+1:3*end/4),yRef(end/2+1:3*end/4) );

xLag = [xL_Ref_1 xR_Ref_1 xL_Ref_2 xR_Ref_2];
yLag = [yL_Ref_1 yR_Ref_1 yL_Ref_2 yR_Ref_2];

plot(xLag,yLag,'ro'); hold on;
axis([0 1 0 1]);

% Combine into ONE Vector
%xLag = [xLag_L xLag_R];
%yLag = [yLag yLag];

% Plot the Geometry
% plot(xLag,yLag,'*'); hold on;
% axis([0 L 0 L]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: rotate geometry
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,y] = rotate_Geometry(ang,xC,yC,xRef,yRef)

len = length(xRef);
x = zeros(1,len); y=x;

xRef = xRef - xC;
yRef = yRef - yC;

for i=1:len
   x(i) =  xRef(i)*cos(ang) - yRef(i)*sin(ang);
   y(i) =  xRef(i)*sin(ang) + yRef(i)*cos(ang);
end

x = x + xC;
y = y + yC;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in the vertices for one coral polyp
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,y,N] = please_Give_Single_Polyp_Geometry(ratio)

xoffset = 0.3; % To put geometry into QUADRANT-1
yoffset = 0.5; % To put geometry into QUADRANT-1

% Note: (1) ds for stored data is 0.6/(2*1024)
%       (2) 'ratio' is comparing 1024:desired resolution

% Get LEFT side geometry
struct_name1 = 'coral2d_left_1024';
[~,x1,y1] = read_Vertex_Points(struct_name1);

% Reverse order so first point is on bottom
x1 = x1(end:-1:1); y1 = y1(end:-1:1);

% Get RIGHT side geometry
struct_name2 = 'coral2d_right_1024';
[~,x2,y2] = read_Vertex_Points(struct_name2);

% Reverse order so first point is on bottom
x2 = x2(end:-1:1); y2 = y2(end:-1:1);

% Put Geometry Together for One Polyp
xAux = [x1; x2]; yAux = [y1; y2];
xAux = xAux+xoffset; yAux = yAux+yoffset;

% Pull out correct resolution
x=xAux(1:ratio:end)';
y=yAux(1:ratio:end)';
N = length(x)/2;


%plot(xAux,yAux,'*'); hold on;
%plot(x,y,'ro'); hold on;

%NOTE: (1) N here is the # of pts. on one arm (NOT ENTIRE POLYP)!!
%      (2) The geometry here is listed for 1024x1024 meshes -> will need
%          to take every 8 or 16 pts to render geometry usable for MATLAB
%          code


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


