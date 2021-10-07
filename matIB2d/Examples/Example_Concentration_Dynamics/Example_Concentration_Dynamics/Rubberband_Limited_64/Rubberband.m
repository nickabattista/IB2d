%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  battistn@tcnj.edu
% Date Modified: April 2021
% Current Institution: TCNJ
% Date Created: May 27th, 2015
% Institution Created: UNC-CH
%
% This code is capable of creating Lagrangian Structures using:
% 	1. Springs
% 	2. Beams (torsional springs or non-invariant beams)
% 	3. Target Points
%	4. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-Tension)")
%
% One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting lengths, etc
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
% 
% If you would like us to add a specific muscle model, please let Nick (battistn@tcnj.edu) know.
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates the RUBBERBAND-EXAMPLE geometry and prints associated input files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Rubberband(Nx,Ny)

%------------------------------------------------------------
% Grid Parameters (MAKE SURE MATCHES IN input2d !!!)
%------------------------------------------------------------
%Nx =  128;      % # of Eulerian Grid Pts. in x-Direction (MUST BE EVEN!!!)
%Ny =  128;      % # of Eulerian Grid Pts. in y-Direction (MUST BE EVEN!!!)
Lx = 1.0;        % Length of Eulerian Grid in x-Direction
Ly = 1.0;        % Length of Eulerian Grid in y-Direction
dx=Lx/Nx;


%---------------------------------------------------------------
% Necessary grid info for ** Concentration ** input file
%  (doesn't need to be changed unless rectangular grid is used)
%---------------------------------------------------------------
x=0:dx:Lx-dx;
y=x;

%------------------------------------------------------------
% Immersed Structure Geometric / Dynamic Parameters %
%------------------------------------------------------------
a = 0.4;         % Length of semi-major axis.
b = 0.2;         % Length of semi-minor axis.
ds_Rest = 0;     % Resting length of springs

%------------------------------------------------------------
% Rubberband overall perimeter length and Lagrangian spacing
%------------------------------------------------------------
L=1.93768964;
ds=Lx/(2*Nx);
s = 0:ds:L; 
N=length(s);

%------------------------------------------------------------
% Structure name
%------------------------------------------------------------
struct_name = 'rubberband'; % Name for .vertex, .spring, etc files.


%------------------------------------------------------------
% Call function to construct geometry
%------------------------------------------------------------
[xLag,yLag] = give_Me_Immsersed_Boundary_Geometry(N,a,b);


%------------------------------------------------------------
% Plot Geometry to test
%------------------------------------------------------------
plot(xLag,yLag,'r-'); hold on;
plot(xLag,yLag,'*'); hold on;
xlabel('x'); ylabel('y');
axis square;

%------------------------------------------------------------
% Prints .vertex file!
%------------------------------------------------------------
print_Lagrangian_Vertices(xLag,yLag,struct_name);


%------------------------------------------------
% Prints .spring file!
%       ->> find k_Spring for certain resolution
%------------------------------------------------
ds128 = 0.5*Lx/128;   % Lag. Spacing in 128x128 Resolution
F_Spr_128 = 4e5;      % Force used in 128x128 Resolution
%
% Get Spring Stiffness Coefficient void of any particular resolution
k_Spring = F_Spr_128 * ds128^3 / ds128;
%
% Get Spring Force based on particular resolution
F_Spring = k_Spring * ds / ds^3;
%
print_Lagrangian_Springs(xLag,yLag,F_Spring,ds_Rest,struct_name);


%------------------------------------------------------------
% Call function for initial concentration
%------------------------------------------------------------
[Concentration,X,Y]= give_Me_Initial_Concentration(Lx,Ly,Nx,Ny,dx,x,y);

%------------------------------------------------------------
% Prints .concentration file!
%------------------------------------------------------------
print_Concentration_Info(Nx,Nx,Concentration,struct_name);

%------------------------------------------------------------
% Prints .geo_connect file!
%------------------------------------------------------------
print_Geometry_Connections(xLag,struct_name);


% Prints .beam file!
%k_Beam = 0.5; C = 0.0;
%print_Lagrangian_Beams(xLag,yLag,k_Beam,C,struct_name);


% Prints .target file!
%k_Target = 1e7;
%print_Lagrangian_Target_Pts(xLag,k_Target,struct_name);

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
% FUNCTION: prints GEOMETRY CONNECTIONS to a file called 
%                               'struct_name'.geo_connect
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Geometry_Connections(xLag,struct_name)


    N = length(xLag); % Total # of Lag. Pts

    geo_fid = fopen([struct_name '.geo_connect'], 'w');

    % Loops over all Lagrangian Pts.
    for i = 1:N
        
        if i<N
            
            s1 = i;   
            s2 = i+1;
            fprintf(geo_fid, '%d %d\n', s1, s2);
            fprintf(geo_fid, '%d %d\n', s2, s1);
            
        elseif i==N
            
            s1 = N;
            s2 = 1;
            fprintf(geo_fid, '%d %d\n', s1, s2);
            fprintf(geo_fid, '%d %d\n', s2, s1);
            
        end
        
    end

    fclose(geo_fid);
        
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints TARGET points to a file called rubberband.vertex
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
% FUNCTION: prints BEAM (Torsional Spring) points to a file called rubberband.beam
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
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',s-1, s, s+1, k_Beam, C);  
            else
                %Case s=N
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',s-1, s, 1,   k_Beam, C);  
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

    fprintf(spring_fid, '%d\n', N );

    %spring_force = kappa_spring*ds/(ds^2);

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

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints initial CONCENTRATION to rubberband.concentration
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
% FUNCTION: creates the Lagrangian structure geometry
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag] = give_Me_Immsersed_Boundary_Geometry(N,rmin,rmax)

% The immsersed structure is an ellipse %
for i=1:N
    
    xLag(i) = 0.5 + rmax * cos( 2*pi/N*(i-1) );
    yLag(i) = 0.5 + rmin * sin( 2*pi/N*(i-1) );
    
end


    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives initial concentration in Eulerian domain
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C,X,Y] = give_Me_Initial_Concentration(Lx,Ly,Nx,Ny,dx,x,y)

[X,Y]=meshgrid(x,y);

C=zeros(Nx,Nx);


