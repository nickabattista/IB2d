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

function VIV_Cylinder2_Geometry()

%
% Grid Parameters (MAKE SURE MATCHES IN input2d !!!)
%
Nx =  256;        % # of Eulerian Grid Pts. in x-Direction (MUST BE EVEN!!!)
Ny =  256;        % # of Eulerian Grid Pts. in y-Direction (MUST BE EVEN!!!)
Lx = 1.0;        % Length of Eulerian Grid in x-Direction
Ly = 1.0;        % Length of Eulerian Grid in y-Direction


% Immersed Structure Geometric / Dynamic Parameters %
ds= min(Lx/(2.1*Nx),Ly/(2.1*Ny));  % Lagrangian spacing
L = 0.9*Lx;                    % Length of Channel
w = 0.2*Ly;                    % Width of Channel
x1_0 = 0.3;                    % x-Center for Cylinder 1
y1_0 = 0.5;                    % y-Center for Cylinder 1
r1 = w/6;                      % Radii of Cylinder 1
x2_0 = 0.5;
y2_0 = 0.5;
r2 = w/6;

struct_name = 'viv_geo2cyl';       % Name for .vertex, .spring, etc files.


% Call function to construct geometry
[xLag,yLag] = give_Me_Channel_Immsersed_Boundary_Geometry(ds,L,w,Lx,Ly);
[x1Lag_Cy,y1Lag_Cy] = give_Me_Cylinder_Immsersed_Boundary_Geometry(ds,r1,x1_0,y1_0);
[x2Lag_Cy,y2Lag_Cy] = give_Me_Cylinder_Immsersed_Boundary_Geometry(ds,r2,x2_0,y2_0);

x1Tether = x1_0;   % xValues you want to tether to on the channel
x2Tether = x2_0;
[indsTether_CY1, x1_0_new] = give_Me_Tethering_Pt_Indices(xLag,x1Tether);
[indsTether_CY2, x2_0_new] = give_Me_Tethering_Pt_Indices(xLag, x2Tether);
x1Lag_Cy = x1Lag_Cy + (x1_0_new-x1_0); % shift circle over to match indices better
x2Lag_Cy = x2Lag_Cy + (x2_0_new - x2_0);

% Plot Geometry to test
plot(xLag(1:end/2),yLag(1:end/2),'r-'); hold on;
plot(xLag(end/2+1:end),yLag(end/2+1:end),'r-'); hold on;

N = length(x1Lag_Cy);
M = length(x2Lag_Cy);

plot(x1Lag_Cy,y1Lag_Cy,'r*'); hold on;
plot(x2Lag_Cy, y2Lag_Cy, 'y*'); hold on;
plot(x1Lag_Cy(1),y1Lag_Cy(1),'g*'); hold on;
plot(x2Lag_Cy(1), y2Lag_Cy(2), 'b*'); hold on;
plot(xLag(indsTether_CY1(1)),yLag(indsTether_CY1(1)),'b*'); hold on;
plot(xLag(indsTether_CY2(1)), yLag(indsTether_CY2(1)), 'r*'); hold on;
plot(xLag(indsTether_CY1(2)),yLag(indsTether_CY1(2)),'g*'); hold on;
plot(xLag(indsTether_CY2(2)), yLag(indsTether_CY2(2)), 'm*'); hold on;
plot(x1Lag_Cy(N/2+1),y1Lag_Cy(N/2+1),'b*'); hold on;
plot(x2Lag_Cy(M/2+1), y2Lag_Cy(M/2+1), 'r*'); hold on;
xlabel('x'); ylabel('y');
axis square;


% Prints .vertex file!
x_bothLag_Cy = [x1Lag_Cy x2Lag_Cy]; % concatenates the two cylinder x coordinate arrays
y_bothLag_Cy = [y1Lag_Cy y2Lag_Cy];
print_Lagrangian_Vertices([xLag x_bothLag_Cy],[yLag y_bothLag_Cy],struct_name, 'w');

% Prints .spring file!
k_Spring = 2.0e7; 
resting_length_tether1 = 2*r1;
offset = length(xLag);
% since cylinders are identical, modified print_Lagrangian_Springs to take
% in data for one, then it prints out two times as much data
print_Lagrangian_Springs(x1Lag_Cy,y1Lag_Cy, k_Spring,ds,r1,offset,indsTether_CY1, resting_length_tether1,struct_name, 'w');


% modified print_Lagrangian_Damped_Springs to take in second tether location
% these are the tethers
% Prints .d_spring file! (FOR DAMPED SPRINGS)
k_Spring = 2e4; 
resting_length_tether1 = 2*r1;
offset = length(xLag); 
b_damp = 5.0;
print_Lagrangian_Damped_Springs(x1Lag_Cy,y1Lag_Cy, k_Spring,ds,r1,offset,indsTether_CY1, indsTether_CY2, resting_length_tether1,b_damp,struct_name, 'w');


% Change this like the undamped spring file
% Prints .beam file!
k_Beam = 5.0e9;  
C1 = compute_Curvatures(x1Lag_Cy,y1Lag_Cy);
print_Lagrangian_Beams(x1Lag_Cy,y1Lag_Cy, k_Beam,C1,struct_name,offset, 'w');

% Prints .target file!
k_Target = 2.5e7;
print_Lagrangian_Target_Pts(xLag,k_Target,struct_name, 'w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints VERTEX points to a file called .vertex
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Vertices(xLag,yLag,struct_name, mode)

    N = length(xLag);

    vertex_fid = fopen([struct_name '.vertex'], mode);

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

function print_Lagrangian_Target_Pts(xLag,k_Target,struct_name, mode)

    N = length(xLag);

    target_fid = fopen([struct_name '.target'], mode);

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

function print_Lagrangian_Beams(xLag,yLag,k_Beam,C,struct_name,offset, mode)

    % k_Beam: beam stiffness
    % C: beam curvature
    
    N = length(xLag); % NOTE: Total number of beams = Number of Total Lag Pts. on Cylinder

    beam_fid = fopen([struct_name '.beam'], mode);

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
    
    %do it again for the second cylinder
    for s = 1:N
            if ( (s>1) && (s <= N-1) )        
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',s-1+offset+N, s+offset+N, s+1+offset+N, k_Beam, C(s) );  
            elseif (s==1)
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',N+N+offset,   N+s+offset, N+s+1+offset, k_Beam, C(s) );  
            elseif (s==N)
                %Case s=N
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',N+N-1+offset, N+N+offset, N+1+offset,   k_Beam, C(s) );  
            end
    end
    
    
    
    fclose(beam_fid); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints SPRING points to a file called .spring
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Springs(xLag,yLag,k_Spring,ds_Rest,r,offset,indsTether,resting_length_tether,struct_name, mode)

    % offset = # of lag pts in channel (ordering lag pts, channel -> cylinder)

    N = length(xLag); % Number of pts on cylinder

    spring_fid = fopen([struct_name '.spring'], mode);

    fprintf(spring_fid, '%d\n', (N+N/2)*2 ); % + 2 if non-damped springs connecting cylinder to channel

    %spring_force = kappa_spring*ds/(ds^2);

    %SPRINGS BETWEEN 1st CYL VERTICES
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
    
    % SPRINGS BETWEEN 2nd CYL VERTICES
    for s = 1:N
            if s < N    
                ds_Rest = sqrt( ( xLag(s) - xLag(s+1) )^2 +  ( yLag(s) - yLag(s+1) )^2 ); 
                fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s+offset+N, s+1+offset+N, k_Spring, ds_Rest);  
            else
                %Case s=N
                ds_Rest = sqrt( ( xLag(s) - xLag(1) )^2 +  ( yLag(s) - yLag(1) )^2 ); 
                fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s+offset+N, 1+offset+N,   k_Spring, ds_Rest);  
            end
    end
    
    for s=1:N/2
        ds_Rest = sqrt( ( xLag(s) - xLag(s+N/2) )^2 +  ( yLag(s) - yLag(s+N/2) )^2 ); 
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s+N+offset, s+N+N/2+offset, k_Spring, ds_Rest);
    end
    
    % UNCOMMENT IF UNDAMPED SPRINGS CONNECT CYLINDER TO CHANNEL
    %s=1; % Reset
    %fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s+offset,     indsTether(2), 1e4, resting_length_tether);
    %fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s+N/2+1+offset, indsTether(1), 1e4, resting_length_tether);
    
    fclose(spring_fid); 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints DAMPED SPRING points to a file called struct.d_spring
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Damped_Springs(xLag,yLag,k_Spring,ds_Rest,r,offset,indsTether1, indsTether2,resting_length_tether,b_damp,struct_name, mode)
    % how function is called above:
    % print_Lagrangian_Damped_Springs(x1Lag_Cy,y1Lag_Cy, k_Spring,ds,r1,
    % offset,indsTether_CY1, resting_length_tether1,b_damp,struct_name, 'w');
    

    N = length(xLag); % just two damped springs to connect cylinder to channel

    dSpring_fid = fopen([struct_name '.d_spring'], mode);

    fprintf(dSpring_fid, '%d\n', 2 );
    
    s=1; % Reset
    fprintf(dSpring_fid, '%d %d %1.16e %1.16e %1.16e\n', s+offset,       indsTether1(2), k_Spring, resting_length_tether,b_damp);
    fprintf(dSpring_fid, '%d %d %1.16e %1.16e %1.16e\n', s+N/2+1+offset, indsTether1(1), k_Spring, resting_length_tether,b_damp);
    
    % DO it again for the second cylinder
    fprintf(dSpring_fid, '%d %d %1.16e %1.16e %1.16e\n', s+offset+N,       indsTether2(2), k_Spring, resting_length_tether,b_damp);
    fprintf(dSpring_fid, '%d %d %1.16e %1.16e %1.16e\n', s+N+N/2+1+offset, indsTether2(1), k_Spring, resting_length_tether,b_damp);
    
    fclose(dSpring_fid);

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

   
