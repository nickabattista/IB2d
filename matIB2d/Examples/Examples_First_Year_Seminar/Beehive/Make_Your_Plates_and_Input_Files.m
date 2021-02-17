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
struct_name = 'hive';     % Name for .vertex, .spring, etc files. (must match what's in 'input2d')


% Call function to construct geometry
[xLag,yLag,dist] = give_Me_Immsersed_Boundary_Geometry(ds,Nx,Lx);


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

% Prints .concentration file!
C_Left = give_Me_Initial_Concentration(Lx,Lx,Nx,Nx,dx,dx,'left');
C_Right = give_Me_Initial_Concentration(Lx,Lx,Nx,Nx,dx,dx,'right');
Concentration = C_Left+C_Right;
print_Concentration_Info(Nx,Nx,Concentration,struct_name);

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
% FUNCTION: prints CONCENTRATION INFO to file called
%           'struct_name'.concentration
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function print_Concentration_Info(Nx,Ny,C,struct_name)

    con_fid = fopen([struct_name '.concentration'], 'w');


    for i=1:Ny
        for j=1:Nx
            fprintf(con_fid, '%1.16e ', C(i,j) );
        end
        fprintf(con_fid,'\n');
    end
    
    
    
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
 
function [xLag,yLag,dist] = give_Me_Immsersed_Boundary_Geometry(ds,Nx,L)

% ds: Lagrangian pt. spacing
% Nx: Eulerian grid resolution
% L:  Length of computational domain

yLag = 0:ds:L/5;              % Make plate of length L/10
xLag = zeros(1,length(yLag)); % Make corresponding xPts for Vertical Plate

% Translate points
yLag = yLag + (L/2-L/10);  % Translate them symmetrically

% Make sides
xLag_Left =  xLag + (L/2-L/20);
xLag_Right = xLag + (L/2+L/20);

% Combine Pts Into ONE Vector
xLag = [xLag_Left xLag_Right];
yLag = [yLag yLag];

% Distance between plates
dist = (L/2+L/20) - (L/2-L/20);

% Make a 2ND Set of Plates
xLag_L = xLag - 2.275*dist;
xLag_R = xLag + 2.275*dist;

% Combine into ONE Vector
xLag = [xLag_L xLag_R];
yLag = [yLag yLag];

% Plot the Geometry
% plot(xLag,yLag,'*'); hold on;
% axis([0 L 0 L]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives initial concentration gradient inside channel
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C = give_Me_Initial_Concentration(Lx,Ly,Nx,Ny,dx,dy,strFlag)

%WHERE OUTER TUBE LIES
%xMin = 0.15; xMax = 0.45;
%yMin = 0.85; yMax = 1.15;

if strcmp(strFlag,'left')
    xMin = 0.225; xMax = 0.32;
    yMin = 0.41; yMax = 0.59;
else % RIGHT SIDE
    xMin = 0.68; xMax = 0.775;
    yMin = 0.41; yMax = 0.59;
end

xMid = (xMin+xMax)/2;
yMid = (yMin+yMax)/2;

xDiff = (xMax-xMin)/2;
yDiff = (yMax-yMin)/2;

x = 0:dx:Lx;
y = 0:dy:Ly;
inds = give_Me_Indices_To_Apply_Force(x,y,xMin,xMax,yMin,yMax);

C = zeros(Ny,Nx);
for i=1:length( inds(:,1) )
    xInd = inds(i,1);
    yInd = inds(i,2);
    xPt = x(xInd);
    yPt = y(yInd);
    %C(xInd,yInd ) = (-0.5/yDiff^2)*( (yPt-yMid) - yDiff )*( (yPt-yMid) + yDiff ) +  (-0.5/xDiff^2)*( (xPt-xMid) - xDiff )*( (xPt-xMid) + xDiff ); %1.0;
    C(yInd,xInd ) = (-1.0/xDiff^2)*( (xPt-xMid) - xDiff )*( (xPt-xMid) + xDiff ); %1.0;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes indices for placing initial concentration 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function inds = give_Me_Indices_To_Apply_Force(x,y,xMin,xMax,yMin,yMax)

j=1; noMinYet = 1;
while noMinYet
    
    if ( x(j) >= xMin )
        iX_min = j;
        noMinYet = 0;
    end
    j=j+1;
end

j=length(x); noMaxYet = 1;
while noMaxYet
    
    if ( x(j) <= xMax )
        iX_max = j;
        noMaxYet = 0;
    end
    j=j-1;
end

j=1; noMinYet = 1;
while noMinYet
    
    if ( y(j) >= yMin )
        iY_min = j;
        noMinYet = 0;
    end
    j=j+1;
end

j=length(y); noMaxYet = 1;
while noMaxYet
    
    if ( y(j) <= yMax )
        iY_max = j;
        noMaxYet = 0;
    end
    j=j-1;
end

iX_Vec = iX_min:1:iX_max;
iY_Vec = iY_min:1:iY_max;

n = 1;
for i=1:length(iX_Vec)
    for j=1:length(iY_Vec)
        inds(n,1) = iX_Vec(i);
        inds(n,2) = iY_Vec(j);
        n = n+1; 
    end
end


