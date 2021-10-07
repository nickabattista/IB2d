%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  battistn@tcnj.edu
% Current Institution: TCNJ
% Date Modified: April 2021
% 
% Date Created: May 27th, 2015
% Institution when Created: UNC-CH
%
% This code is capable of creating Lagrangian Structures using:
% 	1. Springs
% 	2. Beams (*torsional springs or non-invariant beams*)
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
% FUNCTION: creates the LEAF_LEAF-EXAMPLE geometry and prints associated input files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Leaf_Leaf()

%-----------------------------------------------------
% Grid Parameters (MAKE SURE MATCHES IN input2d !!!)
%-----------------------------------------------------
Nx = 384         % # of Eulerian Grid Pts. in x-Direction (MUST BE EVEN!!!)
Ny = 64          % # of Eulerian Grid Pts. in y-Direction (MUST BE EVEN!!!)
Lx = 30.0        % Length of Eulerian Grid in x-Direction (cm)
Ly = 5.0         % Length of Eulerian Grid in y-Direction (cm)
ds= 0.5*Lx/Nx;   % Lagrangian spacing


%----------------------------------------------------------------
% AIR PROPERTIES FOR INPUT
%
%     Density:     [rho] = kg/m^3
%     Newton:      [N] = kg * m/s^2
%     Dyn. Visc.:  [N*s/m^2] = (kg*m/s^2)*(s/m^2) = kg / (m*s)
%       
%----------------------------------------------------------------
% PROPERTIES IN MKS UNITS
%mu_Air = 18.03e-6; % N*s/m^2 dynamic viscosity of air at 18C (64.4F)
%rho_Air = 1.21;    % kg/m^3  density of air at 18C
mu_Air = 1.895e-5; % N*s/m^2 =  dynamic viscosity of air at 18C (64.4F)
rho_Air = 1.145;    % kg/m^3  density of air at 35C
%(from https://www.engineersedge.com/physics/viscosity_of_air_dynamic_and_kinematic_14483.htm)
%
% CONVERT TO CGS UNITS
mu_Air = mu_Air * (1000/1) * (1/100)     %(1000g/1kg) & (1m/100cm) 
rho_Air = rho_Air * (1000/1) * (1/100)^3 %(1000g/1kg) & (1m/100cm)  
%
% FINAL UNITS
% [mu_Air] = g/(cm*s)
% [rho_Air] = g/cm^3


%----------------------------------------------------
% LEAF GEOMETRY
%----------------------------------------------------
xLag = [];                     % Initialize Storage
yLag = [];                     % Initialize Storage
L_leaf = 1.0;                  % Length of tail off cylinder
struct_name = 'leaf'; % Name for .vertex, .spring, etc files
%
x0 = 3.0;                      % Starting x-Pt for LEAF
y0 = Ly/2;                     % Starting y_Pt for LEAF
%
xLag_L = x0:ds:x0+L_leaf;
yLag_L = y0*ones(1,length(xLag_L));


%----------------------------------------------------
% IF CHANNEL WALLS / CYLINDRICAL LEADING EDGE %
%----------------------------------------------------
r = 0.01;                      % Radii of Cylinder

%--------------------------------------------
% Call function to construct geometry
%       * No Channel Walls
%       * No Cylinder on Leading Edge
%--------------------------------------------
[xLag_Cy,yLag_Cy] = give_Me_Cylinder_Immersed_Boundary_Geometry(ds,r,x0,y0);
% 
% % Plot Geometry to test
% plot(xLag_Ch(1:end/2),yLag_Ch(1:end/2),'r-'); hold on;
% plot(xLag_Ch(end/2+1:end),yLag_Ch(end/2+1:end),'r-'); hold on;
% plot(xLag_Cy,yLag_Cy,'r-'); hold on;
% plot(xLag_L,yLag_L,'m-'); hold on;
% axis([0 Lx 0 Ly]);
% %
% plot(xLag_Ch,yLag_Ch,'*'); hold on;
% plot(xLag_Cy,yLag_Cy,'g*'); hold on;
% xlabel('x'); ylabel('y');
% axis square;
% 
% xLag = [xLag xLag_Ch xLag_Cy]; % Add xLagPts from Circle to xLag Pt. Vector (*no springs or beams*)
% yLag = [yLag yLag_Ch yLag_Cy]; % Add xLagPts from Circle to xLag Pt. Vector (*no springs or beams*)
% 
% Nbefore = length(xLag);   % # of Lagrangian Pts Before Leaf to COUNT for...
%                           % ... Cylinder / Channel Points have TARGET PTS
%

xLag = [xLag xLag_Cy]; % Add xLagPts from Circle to xLag Pt. Vector (*no springs or beams*)
yLag = [yLag yLag_Cy]; % Add xLagPts from Circle to xLag Pt. Vector (*no springs or beams*)

Nbefore = length(xLag); % # of Lag Pts Before Leaf Starts
        
xLag = [xLag xLag_L];   % Update xLag Pts to Include Leaf
yLag = [yLag yLag_L];   % Update xLag Pts to Include Leaf

Ntot = length(xLag);    % Total # of Lag Pts

%------------------------
% PLOT GEOMETRY
%------------------------
plot(xLag,yLag,'.'); hold on;
xlabel('x'); ylabel('y');
axis([0 Lx 0 Ly]);

% Prints .vertex file!
print_Lagrangian_Vertices(xLag,yLag,struct_name);


% Prints .spring file!
k_Spring = 2.5e4; %1e7
print_Lagrangian_Springs(xLag,yLag,k_Spring,struct_name,Nbefore,Ntot);


% Prints .beam file!
k_Beam = 2e8; C = 0.0; %2e12 -> 4e12
print_Lagrangian_Beams(xLag,yLag,k_Beam,C,struct_name,Nbefore,Ntot);


% Prints .target file!
k_Target = 2.5e4;
print_Lagrangian_Target_Pts(xLag,k_Target,struct_name,Nbefore);


% Prints .mass file!
k_Mass = 1e0;         % 'spring' stiffness parameter for tethering
Mass = 0.5e-3;        % "MASS" value for 'ghost' lag movement
print_Lagrangian_Mass_Pts(xLag,k_Mass,Mass,struct_name,Nbefore,Ntot);

% Call function for initial concentration
[Concentration]= give_Me_Initial_Concentration(Nx,Ny);

% Prints .concentration file!
print_Concentration_Info(Nx,Ny,Concentration,struct_name);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints VERTEX points to a file called <struct_name>.vertex
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
% FUNCTION: prints TARGET points to a file called <struct_name>.target
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Target_Pts(xLag,k_Target,struct_name,Nbefore)

    %N = length(xLag);

    target_fid = fopen([struct_name '.target'], 'w');

    fprintf(target_fid, '%d\n', Nbefore+3 );

    %Loops over any Channel/Cylinder Pts AND First 3 Pts Along LEAF
    for s = 1:Nbefore+3
        fprintf(target_fid, '%d %1.16e\n', s, k_Target);
    end

    fclose(target_fid); 
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints BEAMs (Torsional Springs) to a file called <struct_name>.beam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Beams(xLag,yLag,k_Beam,C,struct_name,Nbefore,Ntot)
    
    % k_Beam: beam stiffness
    % C: beam curvature
    
    if Nbefore ~= 0
        N = Ntot-Nbefore-1;
    else
        N = Ntot-Nbefore-2; % NOTE: Total number of beams 
    end
    
    
    
    beam_fid = fopen([struct_name '.beam'], 'w');

    fprintf(beam_fid, '%d\n', N );

    %spring_force = kappa_spring*ds/(ds^2);
    
    %------------------------------------------------
    % BEAM BETWEEN CYLINDER AND LEAF
    %       (including those TETHERED in place by 
    %                                   TARGET PTS)
    %------------------------------------------------
    if Nbefore ~= 0 
        s1 = Nbefore;
        s2 = Nbefore+1;
        s3 = Nbefore+2;
        fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',s1, s2, s3, k_Beam, C);   
    end

    %------------------------------------------------
    % BEAM BETWEEN VERTICES ALONG "LEAF"
    %       (including those TETHERED in place by 
    %                                   TARGET PTS)
    %------------------------------------------------
    for s = Nbefore+2:Ntot-1
        s1 = s-1;
        s2 = s;
        s3 = s+1;
        fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',s1, s2, s3, k_Beam, C);  
    end
    fclose(beam_fid); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints SPRING points to a file called <struct_name>.spring
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Springs(xLag,yLag,k_Spring,struct_name,Nbefore,Ntot)

    if Nbefore ~= 0
        N = Ntot-Nbefore;
    else
        N = Ntot-Nbefore-1;
    end

    spring_fid = fopen([struct_name '.spring'], 'w');

    fprintf(spring_fid, '%d\n', N );

    %spring_force = kappa_spring*ds/(ds^2);
    
    %------------------------------------------------
    % SPRINGS BETWEEN CYLINDER AND LEAF
    %       (including those TETHERED in place by 
    %                                   TARGET PTS)
    %------------------------------------------------
    if Nbefore ~= 0
        
        % springs between adjacent pts on cylinder
        for i=1:Nbefore
            if i<Nbefore
                x1 = xLag(i); y1 = yLag(i);
                x2 = xLag(i+1); y2 = yLag(i+1);
                ds = sqrt( (x1-x2)^2 + (y1-y2)^2 );
                fprintf(spring_fid, '%d %d %1.16e %1.16e\n', i, i+1, k_Spring, ds);    
            else
                x1 = xLag(i); y1 = yLag(i);
                x2 = xLag(1); y2 = yLag(1);
                ds = sqrt( (x1-x2)^2 + (y1-y2)^2 );
                fprintf(spring_fid, '%d %d %1.16e %1.16e\n', i, 1, k_Spring, ds);    
                
            end
        end
        
        % springs across opposite sides of cylinder
        for i=1:floor(Nbefore/2)
            x1 = xLag(i); y1 = yLag(i);
            x2 = xLag( i+floor(Nbefore/2) ); y2 = yLag( i+floor(Nbefore/2) );
            ds = sqrt( (x1-x2)^2 + (y1-y2)^2 );
            fprintf(spring_fid, '%d %d %1.16e %1.16e\n', i, i+floor(Nbefore/2), k_Spring, ds);
        end
        
        % spring between cylinder and leaf
        x1 = xLag(Nbefore); y1 = yLag(Nbefore);
        x2 = xLag(Nbefore+1); y2 = yLag(Nbefore+1);
        ds = sqrt( (x1-x2)^2 + (y1-y2)^2 );
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', Nbefore, Nbefore+1, k_Spring, ds);
        
        
    end
    
    %------------------------------------------------
    % SPRINGS BETWEEN VERTICES ALONG LEAF 
    %       (including those TETHERED in place by 
    %                                   TARGET PTS)
    %------------------------------------------------
    for s = Nbefore+1:Ntot-1
        
        x1 = xLag(s); y1 = yLag(s);
        x2 = xLag(s+1); y2 = yLag(s+1);
        
        ds = sqrt( (x1-x2)^2 + (y1-y2)^2 );
        
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, s+1, k_Spring, ds);
        
    end
    fclose(spring_fid); 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints TARGET points to a file called <struct_name>.mass
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
function print_Lagrangian_Mass_Pts(xLag,kMass,Mass,struct_name,Nbefore,Ntot)

    %LOOP OVER LAG PTS FOR MASS PTS IN LAGRANGIAN INDEXING
    N = (Ntot-Nbefore);
    
    mass_fid = fopen([struct_name '.mass'], 'w');

    fprintf(mass_fid, '%d\n', N );

    %Loops over all Lagrangian Pts.
    for s = Nbefore+1:Ntot%Nbefore:Ntot
        fprintf(mass_fid, '%d %1.16e %1.16e\n', s, kMass, Mass );
    end
 
    fclose(mass_fid); 
      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates the Lagrangian structure geometry for channel
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag] = give_Me_Channel_Immersed_Boundary_Geometry(ds,L,w,Lx,Ly)

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

function [xLagN,yLagN] = give_Me_Cylinder_Immersed_Boundary_Geometry(ds,r,x0,y0)

% The immsersed structure is a cylinder %

dtheta = ds/ (4*r);
theta = 0; i=1;
while theta < 2*pi
   xLag(i) = (x0-r-ds) + r*cos(theta);
   yLag(i) = y0 + r*sin(theta);
   theta = theta + dtheta;
   i=i+1;
end

for i=1:length(xLag)
   xLagN(i) = xLag(end-(i-1));
   yLagN(i) = yLag(end-(i-1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints initial concentration to a file called
% <struct_name>.concentration
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

function [C] = give_Me_Initial_Concentration(Nx,Ny)

C = zeros(Ny,Nx);



   
   
