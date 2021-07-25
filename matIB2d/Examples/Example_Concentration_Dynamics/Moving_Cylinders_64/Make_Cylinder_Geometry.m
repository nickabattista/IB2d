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
% FUNCTION: creates the MOVING CYLINDER geometry and prints associated input files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Make_Cylinder_Geometry()

%---------------------------------------------------------
% Grid Parameters (MAKE SURE MATCHES IN input2d !!!)
%---------------------------------------------------------
Nx =  64;       % # of Eulerian Grid Pts. in x-Direction (MUST BE EVEN!!!)
Ny =  64;       % # of Eulerian Grid Pts. in y-Direction (MUST BE EVEN!!!)
Lx = 1.0;        % Length of Eulerian Grid in x-Direction
Ly = 1.0;        % Length of Eulerian Grid in y-Direction
dx = Lx/Nx;
dy = Ly/Ny; 
x = 0:dx:Lx-dx;
y = 0:dy:Ly-dy;

%---------------------------------------------------------
% Immersed Structure Geometric / Dynamic Parameters %
%---------------------------------------------------------
ds = 0.5*Lx/Nx;  % Lagrangian Point Spacing
r = 0.02;       % Radii of Cylinder (mm)
D = 2*r;         % Diamater of Cylinders
w = 3*D + D;     % Spacing between Cylinders (3 Diameters apart + 2 radii's)
struct_name = 'cylinders'; % Name for .vertex, .spring, etc files.

%---------------------------------------------------------
% Call function to construct geometry
%---------------------------------------------------------
[xC,yC,NCircle] = give_Me_Immersed_Boundary_Geometry(ds,r);

%---------------------------------------------------------
% CONSTRUCT MULTIPLE CIRCLES SPACED APART APPROPRIATELY
%           AND TRANSLATED APPROPRIATELY IN DOMAIN
%----------------------------------------------------------
xLag = [xC xC xC xC xC];
yLag = [yC yC+w yC+2*w yC+3*w yC+4*w] + r + ( Ly-(5*D+4*(3*D) ))/2;


%-----------------------------
% Plot Geometry to test
%-----------------------------
plot(xLag,yLag,'r.'); hold on;
plot(xLag,yLag,'*'); hold on;
xlabel('x'); ylabel('y');
axis([0 Lx 0 Ly]);


%-----------------------------
% Make States A --> B
%-----------------------------
xTransFrac = 0.05;
xA = xLag + xTransFrac*Lx;     yA = yLag;
xB = xLag + (1-xTransFrac)*Lx; yB = yLag;


%-----------------------------
% Update_Springs.m Info
%-----------------------------
Re_Desired = 1.5;               % Desired Reynolds number
rho = 1000;                     % Density (kg/m^3)
%
period_Desired = 0.5;           % Finds appropriated velocity
dist = ( (1-2*xTransFrac) )*Lx; % Horizontal Distance Traveled 
speed = dist / period_Desired;  % Period of Movement Each Direction
%
mu = rho*speed*D / Re_Desired;  % Reynolds number
%
fprintf('\n----------------------------------------------------------\n');
fprintf('    >>>> Info for update_Target_Point_Positions.m <<<<\n\n');
fprintf('   -> Circles go from %.2f -> %.2f\n',xTransFrac*Lx,(1-xTransFrac)*Lx);
fprintf('   -> Desired Avg. Horizontal Speed: %.4f\n',speed);
fprintf('   -> Period = %.4f\n',period_Desired);
fprintf('   -> Re = %.4f\n\n',Re_Desired);
%
fprintf('\n----------------------------------------------------------\n');
fprintf('    >>>> Info for input2d <<<<\n\n');
fprintf('   -> mu = %.4f\n\n',mu);




%-----------------------------
% Plot States to rest
%-----------------------------
figure(2)
plot(xA,yA,'b.'); hold on;
plot(xB,yB,'r.'); hold on;
xlabel('x'); ylabel('y');
axis([0 Lx 0 Ly]);

%-----------------------------
% Print STATES A,B, and C
%-----------------------------
print_States(xA,yA,'State_A');
print_States(xB,yB,'State_B');

%-----------------------------
% Prints .vertex file!
%-----------------------------
print_Lagrangian_Vertices(xA,yA,struct_name);

%-----------------------------
% Prints .spring file!
%-----------------------------
k_Spring = 1e6;
print_Lagrangian_Springs(xLag,yLag,k_Spring,NCircle,struct_name);

%-----------------------------
% Prints .beam file!
%-----------------------------
kBeam = 1e9;
print_Lagrangian_nonInv_Beams(xLag,yLag,kBeam,NCircle,struct_name);

%-----------------------------
% Prints .target file!
%-----------------------------
k_Target = 5e6;
print_Lagrangian_Target_Pts(xLag,k_Target,struct_name);

%CONCENTRATION 
% Call function for initial concentration
[Concentration]= give_Me_Initial_Concentration(Nx,Ny,x,y);

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
% FUNCTION: prints BEAM (Torsional Spring) points to a file called <struct_name>.beam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_nonInv_Beams(xLag,yLag,kBeam,NCircle,struct_name)

    % k_Beam: beam stiffness
    % C: beam curvature

    beam_fid = fopen([struct_name '.nonInv_beam'], 'w');

    fprintf(beam_fid, '%d\n', 5*NCircle );

    %beam_force = kappa_beam*ds/(ds^4);

    %---------------------------------------------------------
    % BEAMS BETWEEN ADJACENT VERTICES ON EACH CYLINDER
    %---------------------------------------------------------
    for j=1:5
        
        for s = 1:NCircle
           
            % Auxiliary index to account for 5 cylinders
            ss = s + (j-1)*NCircle; 
           
            if s==1   % Middle Pt is FIRST Point
                
                id_L = j*NCircle; 
                id_M = ss;
                id_R = ss+1;
                
                Cx = xLag(id_L)+xLag(id_R)-2*xLag(id_M);
                Cy = yLag(id_L)+yLag(id_R)-2*yLag(id_M);  
                
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', id_L, id_M, id_R, kBeam, Cx, Cy);  
                
            elseif s==NCircle % Middle Point is LAST Point
                
                id_L = ss-1; 
                id_M = ss;
                id_R = 1+(j-1)*NCircle;
                
                Cx = xLag(id_L)+xLag(id_R)-2*xLag(id_M);
                Cy = yLag(id_L)+yLag(id_R)-2*yLag(id_M);  
                
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', id_L, id_M, id_R, kBeam, Cx, Cy);  
                
                
            else % ALL interior points along each circle
                
                id_L = ss-1; 
                id_M = ss;
                id_R = ss+1;
                
                Cx = xLag(id_L)+xLag(id_R)-2*xLag(id_M);
                Cy = yLag(id_L)+yLag(id_R)-2*yLag(id_M);  
                
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', id_L, id_M, id_R, kBeam, Cx, Cy);  
                 
            end
        end
    end
    
    fclose(beam_fid); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints SPRING points to a file called <struct_name>.spring
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Springs(xLag,yLag,k_Spring,NCircle,struct_name)


    spring_fid = fopen([struct_name '.spring'], 'w');

    fprintf(spring_fid, '%d\n', 5*NCircle + 5*floor(NCircle/2) );

    %spring_force = kappa_spring*ds/(ds^2);

    %---------------------------------------------------------
    % SPRINGS BETWEEN ADJACENT VERTICES ON EACH CYLINDER
    %---------------------------------------------------------
    for j=1:5
        
        for s = 1:NCircle
           
            % Auxiliary index to account for 5 cylinders
            ss = s + (j-1)*NCircle; 
           
            if s < NCircle   
                
                s1 = ss;
                s2 = ss+1;
                ds = sqrt( ( xLag(s1)-xLag(s2) )^2 + ( yLag(s1)-yLag(s2) )^2  );
                
                fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s1, s2, k_Spring, ds);  
            
            else %Case s = NCircle
                
                s1 = ss;
                s2 = 1+(j-1)*NCircle;
                ds = sqrt( ( xLag(s1)-xLag(s2) )^2 + ( yLag(s1)-yLag(s2) )^2  );
                                 
                fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s1, s2, k_Spring, ds);  
            end
        end
    end
    
    %---------------------------------------------------------
    % SPRINGS BETWEEN OPPOSITE SIDE VERTICES ON EACH CYLINDER
    %---------------------------------------------------------
    space = floor(NCircle/2);
    for j=1:5
        
        for s = 1:NCircle
           
            % Auxiliary index to account for 5 cylinders
            ss = s + (j-1)*NCircle; 
           
            if s <= floor(NCircle/2)
                
                s1 = ss;
                s2 = ss+space;
                ds = sqrt( ( xLag(s1)-xLag(s2) )^2 + ( yLag(s1)-yLag(s2) )^2  );
                
                fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s1, s2, k_Spring, ds);  
            
            end
            
        end
    end  
    
    
    
    fclose(spring_fid); 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints VERTEX points to a file called State_<j>.pts
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_States(xA,yA,struct_name)

    N = length(xA);

    vertex_fid = fopen([struct_name '.pts'], 'w');

    %Loops over all Lagrangian Pts.
    for s = 1:N
        X_v = xA(s);
        Y_v = yA(s);
        fprintf(vertex_fid, '%1.16e %1.16e\n', X_v, Y_v);
    end

    fclose(vertex_fid);  
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates the Lagrangian structure geometry
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag,NCircle] = give_Me_Immersed_Boundary_Geometry(ds,r)

Npts = ceil( 2*pi*r / ds );   % Get total # of Lag Pts around Cylinder

tVec = linspace(0,2*pi,Npts); % Vector of Theta values
tVec = tVec(1:end-1);         % Don't need theta=2*pi value at end (overlaps first point)

% The immsersed structure is a circle %
for i=1:length(tVec)
    
    xLag(i) = r * cos( tVec(i) );
    yLag(i) = r * sin( tVec(i) );
    
end

NCircle = length(xLag); % # of pts in circle

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

function [C] = give_Me_Initial_Concentration(Nx,Ny,x,y)

[xx,yy]=meshgrid(x,y);
xc = 0.3; 
C = exp(-7*(2*(xx - xc)./0.4).^2);