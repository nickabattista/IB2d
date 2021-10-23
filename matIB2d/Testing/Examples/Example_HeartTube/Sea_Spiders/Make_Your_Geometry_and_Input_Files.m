%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  nickabattista@gmail.com
% Date Created: January 8th, 2018
% Institution: UNC-CH
%
% This code is capable of creating Lagrangian Structures using:
% 	1. Springs
% 	2. Beams (torsional springs, non-invariant beams)
% 	3. Target Points
%	4. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-Tension)")
%
% One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting lengths, etc
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
% 
% If you would like us to add a specific muscle model, please let Nick (nickabattista@gmail.com) know.
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates the geometry and prints associated input files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Make_Your_Geometry_and_Input_Files()

%
% Grid Parameters (MAKE SURE MATCHES IN input2d !!!)
%
Nx =  1024;      % # of Eulerian Grid Pts. in x-Direction (MUST BE EVEN!!!)
Lx = 1.0;        % Length of Eulerian Grid in x-Direction
dx = Lx/Nx;      % Grid spatial resolution

Ny = 128;        % # of Eulerian Grid Pts. in y-Direction (MUST BE EVEN!!!)
Ly = 0.125;      % Length of Eulerian Grid in y-Direction

%
% Immersed Structure Geometric / Dynamic Parameters %
%
ds = 0.5*dx;                % Lagrangian Pt. Spacing (2x resolution of Eulerian grid)
struct_name = 'sea_spider'; % Name for .vertex, .spring, etc files. (must match what's in 'input2d')
legD = 0.08;                % leg diameter
gutD = 0.5*0.1;             % gut diameter

% Call function to construct geometry
[xLag,yLag,Ninfo] = give_Me_Immsersed_Boundary_Geometry(ds,Nx,Lx,legD,gutD);

% Call function to construct initial background concentration
[Concentration,~] = give_Me_Initial_Concentration(Lx,Ly,Nx,Ny,dx,legD,gutD);


% Plot Geometry to test
plot(xLag,yLag,'r-'); hold on;
plot(xLag,yLag,'*'); hold on;
xlabel('x'); ylabel('y');
axis([0 Lx 0 Lx]);

% Prints USEFUL INFO FOR UPDATE_SPRINGS to screen:
fprintf('\n\n USEFUL FOR UPDATE_SPRINGS: \n');
fprintf(' gutD: %d  (gut diameter)\n',gutD);
fprintf(' Ninfo(1): %d  (# of pts on top of gut)\n',Ninfo(1));
fprintf(' Ninfo(2): %d  (# of pts in gut)\n',Ninfo(2));
fprintf(' Ninfo(3): %d  (# of pts before bottom of gut)\n\n\n',Ninfo(3));

% Prints .vertex file!
print_Lagrangian_Vertices(xLag,yLag,struct_name);


% Prints .spring file!
k_Spring_Adj = 5e7;                   % Spring stiffness (adjacent lag. pts)
k_Spring_Across = 1e3;                % Spring stiffness (across)
ds_Adj = ds;                          % Spring resting length (adjacent)
ds_Across = legD;                     % Spring resting length (across)
print_Lagrangian_Springs(xLag,yLag,k_Spring_Adj,k_Spring_Across,ds_Adj,ds_Across,struct_name,Ninfo)


% Prints .beam (torsional spring) file !
%k_Beam = 1e8;               % Beam Stiffness (does not need to be equal for all beams)
%C = 0;                      % "Curvature" of initial configuration
%print_Lagrangian_Beams(xLag,yLag,k_Beam,C,struct_name,Ninfo);


% Prints .nonInv_beam (torsional spring) file !
k_Beam = 1e7;               % Beam Stiffness (does not need to be equal for all beams)
print_Lagrangian_nonInv_Beams(xLag,yLag,k_Beam,struct_name,Ninfo)


% Prints .target file!
k_Target = 1e6;
print_Lagrangian_Target_Pts(xLag,k_Target,struct_name,Ninfo);

% Prints .porous file!
alpha = 1e-4; 
print_Lagrangian_Porosity(xLag,alpha,struct_name,Ninfo)

% Prints .concentration file!
print_Concentration_Info(Nx,Ny,Concentration,struct_name);


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
% FUNCTION: prints TARGET points to a file called 'struct_name'.target
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Target_Pts(xLag,k_Target,struct_name,Ninfo)

    % Ninfo(2): # of lag. pts. before leg (outer tube) geometry

    N = length(xLag)-Ninfo(2);

    target_fid = fopen([struct_name '.target'], 'w');

    fprintf(target_fid, '%d\n', N );

    % Hold ends of gut in place
    fprintf(target_fid, '%d %1.16e\n', 1, k_Target);
    fprintf(target_fid, '%d %1.16e\n', Ninfo(1), k_Target);
    fprintf(target_fid, '%d %1.16e\n', Ninfo(1)+1, k_Target);    
    fprintf(target_fid, '%d %1.16e\n', Ninfo(2), k_Target);

    
    %Loops over all Lagrangian Pts.
    for s = Ninfo(2)+1:length(xLag)
        fprintf(target_fid, '%d %1.16e\n', s, k_Target);
    end

    fclose(target_fid);   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints POROSITY points to a file called 'struct_name'.porous
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Porosity(xLag,alpha,struct_name,Ninfo)

    % Ninfo(2): # of lag. pts. before leg (outer tube) geometry

    N = length(xLag)-Ninfo(2);
    
    porous_fid = fopen([struct_name '.porous'], 'w');

    fprintf(porous_fid, '%d\n', N );

    %Loops over all Lagrangian Pts.
    for s = Ninfo(2)+1:length(xLag)
        if s == Ninfo(2)+1
            fprintf(porous_fid, '%d %1.16e %1.16e\n', s, alpha,-2);      
        elseif s == Ninfo(2)+2
            fprintf(porous_fid, '%d %1.16e %1.16e\n', s, alpha,-1);
        elseif s== Ninfo(3)-1
            fprintf(porous_fid, '%d %1.16e %1.16e\n', s, alpha,1);
        elseif s==Ninfo(3)
            fprintf(porous_fid, '%d %1.16e %1.16e\n', s, alpha,2);
        elseif s == Ninfo(3)+1
            fprintf(porous_fid, '%d %1.16e %1.16e\n', s, alpha,-2);      
        elseif s == Ninfo(3)+2
            fprintf(porous_fid, '%d %1.16e %1.16e\n', s, alpha,-1);
        elseif s== length(xLag)-1
            fprintf(porous_fid, '%d %1.16e %1.16e\n', s, alpha,1);
        elseif s== length(xLag)
            fprintf(porous_fid, '%d %1.16e %1.16e\n', s, alpha,2);
        else
            fprintf(porous_fid, '%d %1.16e %1.16e\n', s, alpha,0);
        end
    end

    fclose(porous_fid); 

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
% FUNCTION: prints SPRING points to a file called 'struct_name'.spring
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Springs(xLag,yLag,k_Spring_Adj,k_Spring_Across,ds_Adj,ds_Across,struct_name,Ninfo)

    N = Ninfo(2)-2 + Ninfo(1); % adjacent and then across

    spring_fid = fopen([struct_name '.spring'], 'w');

    fprintf(spring_fid, '%d\n', N );    % Print # of springs 

    %SPRINGS ACROSS GUT
    for s = 1:Ninfo(1)
        sTOP = s;
        sBOT = s+Ninfo(1);
        ds_Across = sqrt( (xLag(sTOP)-xLag(sBOT))^2 + (yLag(sTOP)-yLag(sBOT))^2 );
        if ( (s < 50) || ( s > Ninfo(1)-50) )
            fprintf(spring_fid, '%d %d %1.16e %1.16e\n', sTOP, sBOT, 1e-1, ds_Across);
        else
            fprintf(spring_fid, '%d %d %1.16e %1.16e\n', sTOP, sBOT, k_Spring_Across, ds_Across);
        end
    end
    
    
    %SPRINGS ADJACENT VERTICES (top)
    for s = 1:Ninfo(1)-1
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, s+1, k_Spring_Adj, ds_Adj);  
    end
    
    %SPRINGS ADJACENT VERTICES (bot)
    for s = Ninfo(1)+1:Ninfo(2)-1
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, s+1, k_Spring_Adj, ds_Adj);  
    end
    

    fclose(spring_fid);      

    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints BEAM (Torsional Spring) points to a file called 'struct_name'.beam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Beams(xLag,yLag,k_Beam,C,struct_name,Ninfo)

    % k_Beam: beam stiffness
    % C: beam curvature
    
    N = Ninfo(2)-4; % among adjacent lag points

    beam_fid = fopen([struct_name '.beam'], 'w');

    fprintf(beam_fid, '%d\n', N );

    %spring_force = kappa_spring*ds/(ds^2);

    %BEAMS BETWEEN VERTICES (TOP)
    for s = 2:Ninfo(1)-1
        fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',s-1, s, s+1, k_Beam, C );  
    end
    
    %BEAMS BETWEEN VERTICES (BOT)
    for s = Ninfo(1)+2:Ninfo(2)-1
        fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',s-1, s, s+1, k_Beam, C );  
    end
    
   
    fclose(beam_fid); 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints BEAM (NON-INVARIANT) points to a file called 'struct_name'.nonInv_beam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_nonInv_Beams(xLag,yLag,k_Beam,struct_name,Ninfo)

    % k_Beam: beam stiffness
    % C: beam curvature
    
    N = Ninfo(2)-4; % among adjacent lag points

    beam_fid = fopen([struct_name '.nonInv_beam'], 'w');

    fprintf(beam_fid, '%d\n', N );

    %spring_force = kappa_spring*ds/(ds^2);

    %BEAMS BETWEEN VERTICES (TOP)
    for s = 2:Ninfo(1)-1
        fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n',s-1, s, s+1, k_Beam, 0, 0 );  
    end
    
    %BEAMS BETWEEN VERTICES (BOT)
    for s = Ninfo(1)+2:Ninfo(2)-1
        fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n',s-1, s, s+1, k_Beam, 0, 0 );  
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
 
function [xLag,yLag,Ninfo] = give_Me_Immsersed_Boundary_Geometry(ds,Nx,Lx,legD,gutD)
 
% ds:   Lagrangian pt. spacing
% Nx:   Eulerian grid resolution
% legD: leg (tube) diameter
% gutD: gut diameter

% GUT Useful Points
xTubeHor = 0.2*Lx:ds:0.8*Lx;
yTubeHor = 0.5*(1/8)*Lx*ones(1,length(xTubeHor));

%
% GUT (innter tube) Vertices
%
%yGutTop = yTubeHor+gutD/2;
%yGutTop(1:length(xTubeHor)/4) = -0.45*gutD*sin( pi*(xTubeHor(1:length(xTubeHor)/4) - 0.8)/(0.6/4) ) + yTubeHor(1)+gutD/2;
yGutTop = -0.45*gutD*sin( pi*(xTubeHor - 0.2)/(0.6/4) ) + yTubeHor(1)+gutD/2;
yGutBot =  0.45*gutD*sin( pi*(xTubeHor - 0.2)/(0.6/4) ) + yTubeHor(1)-gutD/2;

for i=1:length(yGutTop)
   
    if yGutTop(i)> yTubeHor(1)+gutD/2
        yGutTop(i) =  yTubeHor(1)+gutD/2;
    end
    
    if yGutBot(i) < yTubeHor(1)-gutD/2
        yGutBot(i) =  yTubeHor(1)-gutD/2;
    end
    
    
    
end

%yGutBot = yTubeHor-gutD/2;
%yGutBot(1:length(xTubeHor)/4) = 0.45*gutD*sin( pi*(xTubeHor(1:length(xTubeHor)/4) - 0.8)/(0.6/4) ) + yTubeHor(1)-gutD/2;

xGut = [xTubeHor xTubeHor];
yGut = [yGutTop yGutBot];


%
% LEG Useful Points
%
xTubeHor = 0.1*Lx:8*ds:0.9*Lx;
yTubeHor = 0.5*(1/8)*Lx*ones(1,length(xTubeHor));

%
% LEG (outer tube) Vertices
%
yLegTop = yTubeHor+legD/2;  %+ 0.0025*sin( pi*(xTubeHor - 0.8)/(0.8/16) );
yLegBot = yTubeHor-legD/2;  %- 0.0025*sin( pi*(xTubeHor - 0.8)/(0.8/16) );

%yTubeLeft = 0.5*(1/8)-legD/2+ds:ds:0.5*(1/8)+legD/2-ds;
%xTubeLeft = 0.2*Lx*ones(1,length(yTubeLeft));

% INCLUDING FOOT (left side)
%xLeg = [xTubeHor xTubeHor xTubeLeft];
%yLeg = [yLegTop  yLegBot  yTubeLeft];

% NOT INCLUDING FOOT (left side)
xLeg = [xTubeHor xTubeHor];
yLeg = [yLegTop  yLegBot ];

%
% Combine Geometries
% 
xLag = [xGut xLeg];
yLag = [yGut yLeg];

%
% Number of points info
%
Ninfo(1) = length(xGut)/2; % # pts along top of gut
Ninfo(2) = length(xGut);   % # pts before start of leg
Ninfo(3) = length(xGut)+length(xLeg)/2; % # of points before bottom of leg
% TESTING (i.e., PLOT IT, yo!)
%plot(xGut,yGut,'r*'); hold on
%plot(xLeg,yLeg,'b*'); hold on;
%axis([0 Lx 0 0.125*Lx]);
%pause();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives initial concentration gradient inside channel
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C,inds] = give_Me_Initial_Concentration(Lx,Ly,Nx,Ny,dx,legD,gutD)

%WHERE OUTER TUBE LIES
%xMin = 0.15; xMax = 0.45;
%yMin = 0.85; yMax = 1.15;

xMin = 0.2; xMax = 0.8;

yMid = Ly/2;
yMax = yMid + 1.02*legD/2;
yMin = yMid - 1.02*legD/2;

x = 0:dx:Lx;
y = 0:dx:Ly;

% FOR BOUSSINESQ
%inds = give_Me_Indices_To_Apply_Force(x,y,xMin,xMax,yMin,yMax);
inds = 0;

% initialize
C = zeros(Ny,Nx);

% put in higher concentration
for i=1:length(y)
    for j=1:length(x)
        yVal = y(i);
        xVal = x(j);
        
        if ( ( xVal >= xMin ) && ( xVal <= xMax) )
           
            yMax = 0.5*(1/8)*Lx + legD/2; %+ 0.0025*sin( pi*(xVal - 0.8)/(0.8/16) );
            yMin = 0.5*(1/8)*Lx - legD/2; %- 0.0025*sin( pi*(xVal - 0.8)/(0.8/16) );

            if yVal >= yMax
                C(i,j) = 1; 
            elseif yVal <= yMin
                C(i,j) = 1;
            end
            
           %if ( ( yVal >= yMax ) || ( yVal <= yMin) )
           %     C(i,j) = 1;
           %end
        end
        
    end
end




