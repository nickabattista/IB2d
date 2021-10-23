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

function VIV_Cylinder_Geometry()

%
% Grid Parameters (MAKE SURE MATCHES IN input2d !!!)
%
Nx =  128;        % # of Eulerian Grid Pts. in x-Direction (MUST BE EVEN!!!)
Ny =  128;        % # of Eulerian Grid Pts. in y-Direction (MUST BE EVEN!!!)
Lx = 1.0;        % Length of Eulerian Grid in x-Direction
Ly = 1.0;        % Length of Eulerian Grid in y-Direction


% Immersed Structure Geometric / Dynamic Parameters %
ds= min(Lx/(2*Nx),Ly/(2*Ny));  % Lagrangian spacing
L = 0.9*Lx;                    % Length of Channel
w = 0.2*Ly;                    % Width of Channel
x0 = 0.3;                      % x-Center for Cylinder
y0 = 0.5;                      % y-Center for Cylinder
r = w/6;                       % Radii of Cylinder
struct_name = 'viv_geo';       % Name for .vertex, .spring, etc files.


% Call function to construct geometry
[xLag,yLag] = give_Me_Channel_Immsersed_Boundary_Geometry(ds,L,w,Lx,Ly);
[xLag_Cy,yLag_Cy] = give_Me_Cylinder_Immsersed_Boundary_Geometry(ds,r,x0,y0);

xTether = x0;   % xValues you want to tether to on the channel
[indsTether,x0_new] = give_Me_Tethering_Pt_Indices(xLag,xTether);
xLag_Cy = xLag_Cy + (x0_new-x0); % shift circle over to match indices better

% Plot Geometry to test
plot(xLag(1:end/2),yLag(1:end/2),'r-'); hold on;
plot(xLag(end/2+1:end),yLag(end/2+1:end),'r-'); hold on;
plot(xLag_Cy,yLag_Cy,'r-'); hold on;

N = length(xLag_Cy);
%yLag_Cy(1)
%yLag(indsTether(1))
%yLag_Cy(N/2+1)
%yLag(indsTether(2))

N = length(xLag_Cy);
plot(xLag_Cy,yLag_Cy,'r*'); hold on;
plot(xLag_Cy(1),yLag_Cy(1),'g*'); hold on;
plot(xLag(indsTether(1)),yLag(indsTether(1)),'b*'); hold on;
plot(xLag(indsTether(2)),yLag(indsTether(2)),'g*'); hold on;
plot(xLag_Cy(N/2+1),yLag_Cy(N/2+1),'b*'); hold on;
%plot(xLag_Cy,yLag_Cy,'g*'); hold on;
xlabel('x'); ylabel('y');
axis square;


% Prints .vertex file!
print_Lagrangian_Vertices([xLag xLag_Cy],[yLag yLag_Cy],struct_name);


% Prints .spring file!
k_Spring_Tether = 1e4;
k_Spring = 2.0e7; resting_length_tether = 2*r;
offset = length(xLag);
print_Lagrangian_Springs(xLag_Cy,yLag_Cy,k_Spring,ds,r,offset,indsTether,resting_length_tether,k_Spring_Tether,struct_name);


% Prints .beam file!
k_Beam = 5.0e9;  
C = compute_Curvatures(xLag_Cy,yLag_Cy);
print_Lagrangian_Beams(xLag_Cy,yLag_Cy,k_Beam,C,struct_name,offset);


% Prints .target file!
k_Target = 2.5e7;
print_Lagrangian_Target_Pts(xLag,k_Target,struct_name);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints VERTEX points to a file called .vertex
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
% FUNCTION: prints TARGET points to a file called .target
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
% FUNCTION: prints BEAM (Torsional Spring) points to a file called <str>.beam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Beams(xLag,yLag,k_Beam,C,struct_name,offset)

    % k_Beam: beam stiffness
    % C: beam curvature
    
    N = length(xLag); % NOTE: Total number of beams = Number of Total Lag Pts. on Cylinder

    beam_fid = fopen([struct_name '.beam'], 'w');

    fprintf(beam_fid, '%d\n', N );

    %spring_force = kappa_spring*ds/(ds^2);

    %BEAMS BETWEEN VERTICES
    for s = 1:N
            if ( (s>1) && (s <= N-1) )        
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',s-1+offset, s+offset, s+1+offset, k_Beam, C(s) );  
            elseif (s==1)
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',N+offset,   s+offset, s+1+offset, k_Beam, C(s) );  
            elseif (s==N)
                %Case s=N
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',N-1+offset, N+offset, 1+offset,   k_Beam, C(s) );  
            end
    end
    fclose(beam_fid); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints SPRING points to a file called .spring
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Springs(xLag,yLag,k_Spring,ds_Rest,r,offset,indsTether,resting_length_tether,k_Spring_Tether,struct_name)

    % offset = # of lag pts in channel (ordering lag pts, channel -> cylinder)

    N = length(xLag); % Number of pts on cylinder

    spring_fid = fopen([struct_name '.spring'], 'w');

    fprintf(spring_fid, '%d\n', N+N/2+2 );

    %spring_force = kappa_spring*ds/(ds^2);

    %SPRINGS BETWEEN VERTICES
    for s = 1:N
            if s < N    
                ds_Rest = sqrt( ( xLag(s) - xLag(s+1) )^2 +  ( yLag(s) - yLag(s+1) )^2 ); 
                fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s+offset, s+1+offset, k_Spring, ds_Rest);  
            else
                %Case s=N
                ds_Rest = sqrt( ( xLag(s) - xLag(1) )^2 +  ( yLag(s) - yLag(1) )^2 ); 
                fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s+offset, 1+offset,   k_Spring, ds_Rest);  
            end
    end
    
    for s=1:N/2
        ds_Rest = sqrt( ( xLag(s) - xLag(s+N/2) )^2 +  ( yLag(s) - yLag(s+N/2) )^2 ); 
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s+offset, s+N/2+offset, k_Spring, ds_Rest);
    end
    
    s=1; % Reset
    fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s+offset,     indsTether(2), k_Spring_Tether, resting_length_tether);
    fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s+N/2+1+offset, indsTether(1), k_Spring_Tether, resting_length_tether);
    
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
theta = -pi/2; i=1;
while theta < 3*pi/2
   xLag(i) = x0 - r*cos(theta);
   yLag(i) = y0 - r*sin(theta);
   theta = theta + dtheta;
   i=i+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives indices of Lagrangian Pt choice for tethering.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
function [inds,x0_new] = give_Me_Tethering_Pt_Indices(xLag,xTether)

i=0; not=1;       %initiate for loop
while not == 1
    i=i+1;
    if xLag(i) > xTether
       inds(1) = i;  % to save index for bottom of channel
       not = 0;      % to leave loop
    end
end

next = length(xLag)/2;  
inds(2) = inds(1)+next; % save index for top of channel
x0_new = xLag(inds(1)); % where center of cylinder should be placed.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes "curvature" of ellipse
% 
% NOTE: not curvature in the traditional geometric sense, in the 'discrete'
% sense through cross product.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C = compute_Curvatures(xLag,yLag)

%a-x component (rmin)
%b-y component (rmax)
%C = ab / ( sqrt( a^2*sin(t)^2 + b^2*cos(t)^2  )  )^3

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

   
