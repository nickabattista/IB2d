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
% FUNCTION: creates the HEART-TUBE-EXAMPLE geometry and prints associated input files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function HeartTube()

%
% Grid Parameters (MAKE SURE MATCHES IN input2d !!!)
%
Nx =  128;        % # of Eulerian Grid Pts. in x-Direction (MUST BE EVEN!!!)
Ny =  128;        % # of Eulerian Grid Pts. in y-Direction (MUST BE EVEN!!!)
Lx = 5.0;        % Length of Eulerian Grid in x-Direction
Ly = 5.0;        % Length of Eulerian Grid in y-Direction


% Immersed Structure Geometric / Dynamic Parameters %
ds = Lx/(2*Nx);  % Lagrangian Spacing
d = 1.0;         % Diameter of the tube
L = 4.0;         % Length of the heart-tube
struct_name = 'HeartTube'; % Name for .vertex, .spring, etc files.


% Call function to construct geometry
[xLag,yLag] = give_Me_Immsersed_Boundary_Geometry(ds,L,d,Lx,Ly);


%
% GIVES INITIAL CONCENTRATION 
%
dx = 2*ds; dy = dx;
Concentration = give_Me_Initial_Concentration(Lx,Ly,Nx,Ny,dx,dy);


% Plot Geometry to test
plot(xLag(1:end/2),yLag(1:end/2),'r-'); hold on;
plot(xLag(end/2+1:end),yLag(end/2+1:end),'r-'); hold on;
plot(xLag,yLag,'*'); hold on;
xlabel('x'); ylabel('y');
axis([0 Lx 0 Ly]);


% Prints .vertex file!
print_Lagrangian_Vertices(xLag,yLag,struct_name);


% Prints .spring file!
k_Spring = 1e6;
print_Lagrangian_Springs(xLag,k_Spring,ds,d,struct_name);

% Prints .muscle file!
%LFO = d; SK = 0.3; a = 0.25; b = 4.0; Fmax = 1e5;
%print_Lagrangian_Muscles(xLag,LFO,SK,a,b,Fmax,struct_name)


% Prints .beam file!
k_Beam = 1e2; C = 0.0;
print_Lagrangian_Beams(xLag,k_Beam,C,struct_name);


% Prints .target file!
k_Target = 1e6;
print_Lagrangian_Target_Pts(xLag,k_Target,struct_name);

% Prints .concentration file!
print_Concentration_Info(Nx,Ny,Concentration,struct_name);

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
% FUNCTION: prints Vertex points to a file called rubberband.vertex
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Target_Pts(xLag,k_Target,struct_name)

    N = length(xLag); % Total Number of Lagrangian Pts
    num = 1;          % Number of target points on each end

    target_fid = fopen([struct_name '.target'], 'w');

    fprintf(target_fid, '%d\n', 4*num );

    %Loops over all Lagrangian Pts.
    %for s = 1:N
    
    
    %Left Bottom Target Points
    for s=1:num
        fprintf(target_fid, '%d %1.16e\n', s, k_Target);
    end
    
    %Right Bottom Target Points
    for s=(N/2)-(num-1):N/2
        fprintf(target_fid, '%d %1.16e\n', s, k_Target);
    end
    
    %Left Top Target Points
    for s=(N/2+1):(N/2+1)+(num-1)
        fprintf(target_fid, '%d %1.16e\n', s, k_Target);
    end
    
    %Right Top Target Points
    for s=N-(num-1):N
        fprintf(target_fid, '%d %1.16e\n', s, k_Target);
    end
        


    %end

    fclose(target_fid); 
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints BEAM (Torsional Spring) points to a file called rubberband.beam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Beams(xLag,k_Beam,C,struct_name)

    % k_Beam: beam stiffness
    % C: beam curvature
    
    N = length(xLag); % NOTE: Total number of beams = Number of Total Lag Pts. - 2

    beam_fid = fopen([struct_name '.beam'], 'w');

    fprintf(beam_fid, '%d\n', N - 4 );

    %spring_force = kappa_spring*ds/(ds^2);

    %BEAMS BETWEEN VERTICES
    for s = 2:N-1
            if  s <= N/2-1
                % Bottom of tube
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',s-1, s, s+1, k_Beam, C);  
            elseif ( ( s >= N/2+2 ) && ( s <= N-1 ) )
                % Top of Tube
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',s-1, s, s+1,   k_Beam, C);  
            end
    end
    fclose(beam_fid); 
    
    
    

    
    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints MUSCLE points to a file called "struct_name".muscle
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Muscles(xLag,LFO,SK,a,b,Fmax,struct_name)

    N = length(xLag); %Number of Lagrangian Pts. Total

    muscle_fid = fopen([struct_name '.muscle'], 'w');

    fprintf(muscle_fid, '%d\n', N/2-2 );

    %spring_force = kappa_spring*ds/(ds^2);

    %MUSCLES BETWEEN VERTICES
    for s = 2:N/2-1
            fprintf(muscle_fid, '%d %d %1.16e %1.16e %1.16e %1.16e %1.16e\n', s, s+N/2, LFO, SK, a,b,Fmax);  
    end
    fclose(muscle_fid);
    
    
    
    
    
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints SPRING points to a file called rubberband.spring
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Springs(xLag,k_Spring,ds_Rest,d,struct_name)

    N = length(xLag); %Number of Lagrangian Pts. Total

    spring_fid = fopen([struct_name '.spring'], 'w');

    fprintf(spring_fid, '%d\n', N-2 + N/2 ); %(N-2) btwn adajcent Lag. Pts, (N/2) btwn opposite sides of HT

    %spring_force = kappa_spring*ds/(ds^2);

    %SPRINGS BETWEEN VERTICES
    for s = 1:N
            if s < N/2         
                fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, s+1, k_Spring, ds_Rest);  
            elseif ( ( s >= N/2+1 ) && ( s < N ) )
                %Case s=N
                fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, s+1, k_Spring, ds_Rest);  
            end
    end
    
    for s = 1:N/2
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, s+N/2,   k_Spring/100, d); 
    end
    fclose(spring_fid); 
    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates the Lagrangian structure geometry
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag] = give_Me_Immsersed_Boundary_Geometry(ds,L,d,Lx,Ly)

% The immsersed structure is a straight heart-tube %
x1 = [-L/2:ds:L/2 L/2];             % Constructs x-Values for bottom of tube
y1 = -d/2*ones(1,length(x1)); % Constructs y-Values for bottom of tube

x2 = x1;                      % Constructs x-Values for top of tube
y2 = -y1;                     % Constructs y-Values for top of tube

x1 = x1 + Lx/2;               % Shift into correct box
x2 = x2 + Lx/2;                % Shift into correct box

y1 = y1 + Ly/2;               % Shift into correct box
y2 = y2 + Ly/2;               % Shift into correct box

xLag = [x1 x2];
yLag = [y1 y2];


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives initial concentration gradient inside channel
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C = give_Me_Initial_Concentration(Lx,Ly,Nx,Ny,dx,dy)

%WHERE OUTER TUBE LIES
%xMin = 0.15; xMax = 0.45;
%yMin = 0.85; yMax = 1.15;

xMin = 1.6; xMax = 2.2;
yMin = 2.1; yMax = 2.9;

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
