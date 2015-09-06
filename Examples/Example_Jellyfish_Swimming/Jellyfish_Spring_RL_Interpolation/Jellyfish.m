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
% FUNCTION: creates the RUBBERBAND-EXAMPLE geometry and prints associated input files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Jellyfish()

%
% Grid Parameters (MAKE SURE MATCHES IN input2d !!!)
%
Nx =  128;        % # of Eulerian Grid Pts. in x-Direction (MUST BE EVEN!!!)
Ny =  128;        % # of Eulerian Grid Pts. in y-Direction (MUST BE EVEN!!!)
Lx = 1.0;        % Length of Eulerian Grid in x-Direction
Ly = 1.0;        % Length of Eulerian Grid in y-Direction


% Immersed Structure Geometric / Dynamic Parameters %
ds = Lx/(2*Nx);  % Lagrangian Pt. Spacing (2x resolution of Eulerian grid)
rmax_1 = 0.5/4;         % Length of semi-major axis for PHASE ONE.
rmin_1 = 0.3/4;         % Length of semi-minor axis for PHASE ONE.
rmax_2 = 0.5/4;         % Length of semi-major axis for PHASE TWO.
rmin_2 = 0.15/4;        % Length of semi-minor axis for PHASE TWO.
struct_name = 'jelly'; % Name for .vertex, .spring, etc files.


% Call function to construct geometry for PHASE ONE
[xLag_1,yLag_1,C1] = give_Me_Immsersed_Boundary_Geometry_P1(ds,rmin_1,rmax_1);
N_lag = length(xLag_1);

% Call function construct geometry for PHASE TWO
[xLag_2,yLag_2,C2,ds2] = give_Me_Immsersed_Boundary_Geometry_P2(N_lag,rmin_2,rmax_2);


% Reflect Geometry
xLag_1 = [xLag_1 -xLag_1(end-1:-1:1)];
yLag_1 = [yLag_1  yLag_1(end-1:-1:1)];
xLag_2 = [xLag_2 -xLag_2(end-1:-1:1)];
yLag_2 = [yLag_2  yLag_2(end-1:-1:1)];


% Translate Geometry
yLag_1 = yLag_1 + 2*rmax_1;
xLag_1 = xLag_1 + Lx/2;
yLag_2 = yLag_2 + 2*rmax_1;
xLag_2 = xLag_2 + Lx/2;

% Gives RESTING-LENGTHS for each phase if interpolating resting lengths for springs
[RL_Bell,RL_Body] = give_Me_Phase_Resting_Lengths(xLag_1,yLag_1,xLag_2,yLag_2);





% Plot Geometry to test
plot(xLag_1,yLag_1,'r-'); hold on;
plot(xLag_1,yLag_1,'*'); hold on;
plot(xLag_2,yLag_2,'r-'); hold on;
plot(xLag_2,yLag_2,'g*'); hold on;
xlabel('x'); ylabel('y');
axis([0 1 0 1]);

% Prints .vertex file!
print_Lagrangian_Vertices(xLag_1,yLag_1,struct_name);


% Prints .spring file!
k_Spring = 1e7;
print_Lagrangian_Springs(xLag_1,yLag_1,k_Spring,struct_name,ds,ds2);


% Prints .beam file!
%k_Beam = 1e9;
%print_Lagrangian_Beams(xLag_1,yLag_1,k_Beam,C1,C2,struct_name);


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
% FUNCTION: prints Vertex points to a file called rubberband.vertex
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

function print_Lagrangian_Beams(xLag,yLag,k_Beam,C,C2,struct_name)

    % k_Beam: beam stiffness
    % C: beam curvature
    
    N = length(xLag); % NOTE: Total number of beams = Number of Total Lag Pts. - 2

    beam_fid = fopen([struct_name '.beam'], 'w');

    fprintf(beam_fid, '%d\n', N-2 );

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

    %fprintf('\nHEADS UP! Print following into the update_Beams file!\n\n');
    %fprintf('\n\nFIRST COL: C1  2ND COL: C2 = \n');
    %for i=1:length(C2)
    %   fprintf('%d %d\n',C(i),C2(i));
    %end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints SPRING points to a file called rubberband.spring
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Springs(xLag,yLag,k_Spring,struct_name,ds1,ds2)

    N = length(xLag);

    N_springs = (N-1);   %(N-1) for around jellybell, (N-1)/2 for btwn bell sides
    
    spring_fid = fopen([struct_name '.spring'], 'w');
    
    fprintf(spring_fid, '%d\n', N_springs );

    %spring_force = kappa_spring*ds/(ds^2);

    %SPRINGS BETWEEN VERTICES
    for s = 1:N
            if s < N         
                x1 = xLag(s);    y1 = yLag(s);
                x2 = xLag(s+1);  y2 = yLag(s+1);
                ds = sqrt( (x1-x2)^2 + (y1-y2)^2 );
                fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, s+1, k_Spring, ds);  
            end
            
            % GEOMETRY IS NOT CLOSED!
            %else
            %    %Case s=N
            %    x1 = xLag(s);   y1 = yLag(s);
            %    x2 = xLag(1);   y2 = yLag(1);
            %    ds = sqrt( (x1-x2)^2 + (y1-y2)^2 );
            %    fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, 1,   k_Spring, ds);  
            %end
            
    end
    
    %SPRINGS BETWEEN ENDS OF JELLYFISH BELL
    for s=1:(N-1)/2
            x1 = xLag(s);    y1 = yLag(s);
            id2 = length(xLag) - (s-1);
            x2 = xLag( id2 );  y2 = yLag( id2 );
            ds(s) = sqrt( (x1-x2)^2 + (y1-y2)^2 );
            if s<=25
                strength = k_Spring/5e4;
            else
                strength = k_Spring/1e2;
            end
            fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, id2, strength, ds(s) );  
    end
    
    
    fprintf('\nHEADS UP! Print following into the update_Springs file!\n\n');
    
    fprintf('N_lagpts = %d (start of springs btwn bell sides)\n\n',N);
    fprintf('ds1 = %d and ds2 = %d\n\n',ds1,ds2);
    fprintf('MAKE SURE TO HARDCODE RESTING LENGTHS FOR BOTH PHASES IN UPDATE!\n\n');
    fprintf('dist_Vector = \n');
    for i=1:25 
       fprintf('%d\n',ds(i));
    end
    fprintf('\n');
    
    fclose(spring_fid); 
    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives actually ELLIPTICAL piece
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,y,t] = compute_ELLIPTIC_Branch(ds,rmin,rmax,ang)

% ang: angle below positive x-axis where ellipse is to start
if ang < 0
    angEnd = abs(ang) + pi;
elseif ang == 0
    angEnd = pi;
end

%initiate
ct = 1;
t(1) = ang;
xN = rmin*cos(ang);   x(ct) = xN;
yN = rmax*sin(ang);   y(ct) = yN;

while t <= pi/2-ds; 
    
    %update counter
    ct = ct+1;
    
    xP = x(ct-1);           % x-Prev
    yP = y(ct-1);           % y-Prev
    
    tN = t(ct-1);                 % previous angle
    tF = tN + pi/20;         % far guess
    tGuess = (tN + tF)/2;   % guess
    
    xN1 = rmin*cos(tGuess);   % x-guess 
    yN1 = rmax*sin(tGuess);   % y-guess
    err = ( ds - sqrt( (xN1-xP)^2 + (yN1-yP)^2 ) );
    
    while abs(err) > 1e-6
       
       if err > 0
          
          tN = tGuess;              % Update 'close' PT. [tN,tGuess,tF]
          tGuess = (tN+tF)/2;       % New Guess
          xN1 = rmin*cos(tGuess);   % x-guess 
          yN1 = rmax*sin(tGuess);   % y-guess          
          
       elseif err < 0
           
          tF = tGuess;             % Update FAR PT. [tN,tGuess,tF] 
          tGuess = (tF+tN)/2;     % New Guess
          xN1 = rmin*cos(tGuess);   % x-guess 
          yN1 = rmax*sin(tGuess);   % y-guess  
          
       end
        
       %compute error
       err = ( ds - sqrt( (xN1-xP)^2 + (yN1-yP)^2 ) );
        
    end

    %save values
    x(ct) = xN1;
    y(ct) = yN1;
    
    %update
    t(ct)=tGuess;
    
end

x(ct) = rmin*cos(pi/2);  
y(ct) = rmax*sin(pi/2);  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates the Lagrangian structure geometry for PHASE ONE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag,C] = give_Me_Immsersed_Boundary_Geometry_P1(ds,rmin,rmax)

% The immsersed structure is HALF an ellipse %

ang = 0;
[xLag,yLag,angs] = compute_ELLIPTIC_Branch(ds,rmin,rmax,ang);

C = compute_Curvatures(xLag,yLag);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: compute ELLIPTIC CIRCUMFERENCE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function arcLength = compute_Ellipse_Circumference(a,b)

h = (a-b)^2 / (a+b)^2;
arcLength = pi*(a+b)*( 1 + 3*h / (10 + sqrt(4-3*h) ) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates the Lagrangian structure geometry for PHASE TWO
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag,C,ds2] = give_Me_Immsersed_Boundary_Geometry_P2(N_lag,rmin2,rmax2)

% The immsersed structure is HALF an ellipse %

arcLength = compute_Ellipse_Circumference(rmin2,rmax2);

ds2 = (arcLength/4) / (N_lag-1);

ang = 0;
[xLag,yLag,angs] = compute_ELLIPTIC_Branch(ds2,rmin2,rmax2,ang);

C = compute_Curvatures(xLag,yLag);

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
for i=2:N-1
   
    %t = angs(i);
   %t = atan( rmin/rmax * tan(t) ); %phi parameter?
   %C(i) = -rmin*rmax / ( sqrt( (rmax*cos(t))^2 + (rmin*sin(t))^2 ) )^3;
   
   % Pts Xp -> Xq -> Xr (same as beam force calc.)
   Xp = xLag(i-1); Xq = xLag(i); Xr = xLag(i+1);
   Yp = yLag(i-1); Yq = yLag(i); Yr = yLag(i+1);
   
   C(i) = (Xr-Xq)*(Yq-Yp) - (Yr-Yq)*(Xq-Xp); %Cross product btwn vectors
   
   
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: finds distances between all two phase for spring resting
%           lengths
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [RL_Bell,RL_Body] = give_Me_Phase_Resting_Lengths(xLag_1,yLag_1,xLag_2,yLag_2)

N = length(xLag_1);
Ntest=length(xLag_2);

if N~=Ntest
   error('NUMBER OF LAG PTS BETWEEN PHASES ARE DIFFERENT!!!!!!');
end

RL_Bell = zeros( (N-1)/2,2);
RL_Body = zeros(N-1,2);

for s=1:(N-1)/2
   
   x1 = xLag_1(s); y1 = yLag_1(s);
   id2 = N - (s-1);
   x2 = xLag_1(id2); y2 = yLag_1(id2);
   dist = sqrt( (x1-x2)^2 + (y1-y2)^2 );
   RL_Bell(s,1) = dist; 
   
   x1 = xLag_2(s); y1 = yLag_2(s);
   id2 = N - (s-1);
   x2 = xLag_2(id2); y2 = yLag_2(id2);
   dist = sqrt( (x1-x2)^2 + (y1-y2)^2 );
   RL_Bell(s,2) = dist; 
   
end

for s=1:N-1
   
   x1 = xLag_1(s); y1 = yLag_1(s);
   x2 = xLag_1(s+1); y2 = yLag_1(s+1);
   dist = sqrt( (x1-x2)^2 + (y1-y2)^2 );
   RL_Body(s,1) = dist; 
   
   x1 = xLag_2(s); y1 = yLag_2(s);
   x2 = xLag_2(s+1); y2 = yLag_2(s+1);
   dist = sqrt( (x1-x2)^2 + (y1-y2)^2 );
   RL_Body(s,2) = dist; 
    
end