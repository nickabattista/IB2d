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
% If you would like us to add a specific muscle model, please let Nick (nick.battista@unc.edu) know.
%
%--------------------------------------------------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates the JELLYFISH-EXAMPLE geometry and prints associated input files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Swimmer_Geometry()


% FLUID GRID PARAMETERS %
L = 8;                              % Length of computational domain (m)
N = 512;                            % # of Cartesian grid meshwidths
dx = L/N;                           % Cartesian mesh width (m)


% Construct Geometry
[xLag,yLag,yLag2,ds] = give_Me_Immsersed_Boundary_Geometry(N,L);

% NAMING CONVENTION FOR SIMULATION 
struct_name = 'swimmer';      % structure name

%
% PRINT INPUT FILES (.vertex, .spring, .beam, etc) %
%

% print vertices
print_Lagrangian_Vertices(xLag,yLag,struct_name);

% print vertex stages for update_curvature
print_Lagrangian_Vertices_Phases(xLag,yLag,yLag2,struct_name);

% print springs
k_Spring = 7.5*1.2750000000000000e+9;        % spring constant (Newton)
print_Lagrangian_Springs(xLag,k_Spring,ds,struct_name);

% print beams
k_Beam = 2.0363359375000002e+12;   % beam stiffness constant (Newton m^2)
print_Lagrangian_nonInv_Beams(xLag,yLag,k_Beam,struct_name);




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
% FUNCTION: prints VERTEX points to a file called rubberband.vertex
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Vertices_Phases(xLag,yLag,yLag2,struct_name)

    N = length(yLag);

    vertex_fid = fopen([struct_name '.phases'], 'w');

    %Loops over all Lagrangian Pts.
    for s = 1:N
        X_v1 = xLag(s);
        Y_v1 = yLag(s);
        Y_v2 = yLag2(s);
        fprintf(vertex_fid, '%1.16e %1.16e %1.16e\n',X_v1,Y_v1, Y_v2);
        %fprintf(vertex_fid,'%1.16e %1.16e\n',Y_v1, Y_v2);
    end

    fclose(vertex_fid);     
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints SPRING points to a file called rubberband.spring
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Springs(xLag,k_Spring,ds_Rest,struct_name)

    N = length(xLag);

    spring_fid = fopen([struct_name '.spring'], 'w');

    fprintf(spring_fid, '%d\n', N-1 );

    %spring_force = kappa_spring*ds/(ds^2);

    % SPRINGS BETWEEN VERTICES ON RHS
    for s = 1:N-1
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, s+1, k_Spring, ds_Rest);  
    end
    
    
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints BEAM (NON-INVARIANT) points to a file called rubberband.nonInv_beam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_nonInv_Beams(xLag,yLag,k_Beam,struct_name)

    % k_Beam: beam stiffness
    % Cx/Cy:  beam curvatures in x/y respestively
    
    N = length(xLag); % NOTE: Total number of beams = Number of Total Lag Pts. - 2

    beam_fid = fopen([struct_name '.nonInv_beam'], 'w');

    fprintf(beam_fid, '%d\n', N-2 );

    % beam_force = kappa_beam*ds/(ds^4)
    
    %BEAMS BETWEEN VERTICES ON RHS
    for s = 1:N-2

        Cx = xLag(s) - 2*xLag(s+1) + xLag(s+2);
        Cy = yLag(s) - 2*yLag(s+1) + yLag(s+2);
        fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', s,s+1,s+2, k_Beam, Cx, Cy);  
 
    end
    

    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates the Lagrangian structure geometry
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag,yLag2,ds] = give_Me_Immsersed_Boundary_Geometry(N,L)

ds = L/(2*N);

% Straight Portion %
xStraight = 0:ds:L/20;
yStraight = zeros(1,length(xStraight));
Nstraight = length(xStraight)-2;


% Curved Portion %

tol = 1e-5;    % error tolerance
n = 1;         % counting for array
xC(n) = 0;     % first x-value at (0,-b)
yC(n) = 0;     % first y-value at (0,-b)
tprev = xC(1); % initial x-value
pow = 3;       % power to create curve

while tprev <= L/10;

    n = n+1;
    tnext = tprev + ds/2; %Guess for next angle value
    tfar = 3*tnext;

    %initiating guess for bisection-algorithm
    xn = tnext;
    yn = tnext^pow;
    errSign = ( ds - sqrt( (xn-xC(n-1))^2 + (yn-yC(n-1))^2 ) );
    err = abs(errSign);

    %Bisection algorithm to make points equally spaced
    while ( err > tol )

        if errSign < 0
            tfar = tnext;
            tnext = (tnext+tprev)/2;
        elseif errSign > 0
            tprev = tnext;
            tnext = (tnext+tfar)/2;
        end

        xn = tnext;
        yn = tnext^pow;
        errSign = ( ds - sqrt( (xn-xC(n-1))^2 + (yn-yC(n-1))^2 ) );
        err = abs(errSign);
    end

    xC(n) = xn;    %Store X-value
    yC(n) = yn;    %Store Y-value
    tprev = tnext; %Update X-value

end

% Save Curved Values % 
xCurve = xC + xStraight(end) + ds;
yCurve = yC;

length(yCurve)

yCurve2 = -yC;

xLag = [xStraight xCurve]; xLag = xLag + 0.1*L; %0.7L
yLag = [yStraight yCurve]; yLag = yLag + 0.5*L;

yLag2 = [yStraight yCurve2]; yLag2 = yLag2 + 0.5*L;

ms = 34;
plot(xLag,yLag,'.','MarkerSize',ms); hold on; 
plot(xLag+L/2,yLag2,'r.','MarkerSize',ms); hold on; 
axis([0 L 0 L]); 

