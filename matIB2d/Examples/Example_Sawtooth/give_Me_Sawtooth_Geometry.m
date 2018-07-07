%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  nickabattista@gmail.com
% Date Created: May 27th, 2015
% Institution: UNC-CH
%
% This code is capable of creating Lagrangian Structures using:
% 	1. Springs
% 	2. Beams (*torsional springs)
% 	3. Target Points
%	4. Muscle-Model (combined Force-Length-Velocity model, "Hill+(Length-Tension)")
%
% One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting lengths, etc
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
% 
% If you would like us to add a specific muscle model, please let Nick (nickabattista@gmail.com) know.
%
%--------------------------------------------------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: updates the target point positions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function give_Me_Sawtooth_Geometry()

Lx = 12;        % Length of Computational Domain in X (m)
Ly = 3.0;        % Length of Computational Domain in Y (m)
Nx = 1024;        % Grid Resolution in X
Ny = 256;        % Grid Resolution in Y
ds = Lx/(2.1*Nx);  % Makes Lagrangian Spaced based on 512x512 finest grid

% Experiment: 
% 12 teeth
% d = 6.3mm (height of teeth)
% L = 7.6cm
% w = 6.3mm
% G = avg. gap between panels
% f = frequency
% a = peak to peak amplitude
%
% RE CALCULATION
% v(characteristic) = afL/G NOTE: a/G = 0cclusion Ratio
% L(characteristic) = d

% Geometric Parameters
d = 0.63; % height of tooth
w = 0.63; % width of tooth
G = 0.325;  % vertical distance between teeth (gap)

%
% Gives Sawtooth Geometry
% 
[xLag, yLag] = give_Me_Lagrangian_Geometry(ds,d,w,G);
xLag = xLag + ( Lx/2- xLag(end/4) ); % shift geometry to center of domain
yLag = yLag + ( Ly/2 );    % shift geometry accordingly to center of domain


%
% PRINTING INFO FOR UPDATE_TARGET_POINT_POSITIONS()
%
fprintf('\n\n Info for Update_Target_Point_Positions():\n\n'); 
fprintf('G = %d (vertical gap btwn top and bottom)\n',G); 
fprintf('N_Top = %d (# of points on top of sawtooth)\n\n', length(xLag)/2);

%
% PLOT FINAL GEOMETRY
% 
ms = 10;
plot(xLag,yLag,'.','MarkerSize',ms); hold on;
axis([0 Lx 0 Ly]);

% Prints .vertex file!
struct_name = 'sawtooth';
print_Lagrangian_Vertices(xLag,yLag,struct_name);

% Prints .target file!
k_Target = 1e11;
print_Lagrangian_Target_Pts(xLag,k_Target,struct_name);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints VERTEX points to a file called sawtooth.vertex
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
% FUNCTION: prints Vertex points to a file called sawtooth.target
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives Lagrangian geometry for sawtooth pump
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[xLag, yLag] = give_Me_Lagrangian_Geometry(ds,d,w,G)

%
% gives geometry for one tooth
%
[xTooth, yTooth] = give_Me_One_Tooth_Geometry(ds,d,w);

%combines for many teeth
numTeeth = 12;
xTLag = xTooth; xBLag = xTooth;
yTLag = yTooth; yBLag = -yTooth;
for j=1:numTeeth-1
    % top portion of sawtooth
    xTLag = [xTLag xTooth+j*w];
    yTLag = [yTLag yTooth];
    
    %bottom portion of sawtooth
    xBLag = [xBLag xTooth+j*w];
    yBLag = [yBLag -yTooth];
end

% testing geometry (when numTeeth=3)
% plot(xLag(1:end/3),yLag(1:end/3),'.'); hold on;
% plot(xLag(end/3+1:2/3*end),yLag(end/3+1:2/3*end),'o'); hold on;

% testing geometry for ALL teeth
%plot(xTLag,yTLag,'r.'); hold on;
%plot(xBLag,yBLag,'b.'); hold on;

%
% COMBINE TOP AND BOTTOM INTO SINGLE VECTOR AND TRANSLATE CORRECT DISTANCE
% (a/2)
%
xLag = [xTLag xBLag];
yLag = [yTLag+G/2 yBLag-G/2];

%Shift
%xLag = xLag + xLag(end/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives Lagrangian geometry for ONE sawtooth 'tooth' geometry
%           NOTE: geometry lower triangular pt: (0,0)
%                 right most point: (w,0)
%                 top-left most point: (0,d)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xTooth, yTooth] = give_Me_One_Tooth_Geometry(ds,d,w)

% Make Vertical Piece
yH = 0:ds:d;
xH = zeros(1,length(yH));

% Make Diagonal Piece
diag = sqrt( d^2 + w^2 );
xDiag = -diag/2:ds:diag/2;
yDiag = zeros(1 , length(xDiag) );

% Rotate Diagonal Piece
ang = atan(d/w);
[xDiag,yDiag] = rotate_geometry(-ang,xDiag,yDiag);

% Translate into correct place
xDiag = xDiag + d/2;
yDiag = yDiag + w/2;

% Shift last point slightly so no overlap when combining teeth
xDiag(end) = ( 0.4*xDiag(end-1)+0.6*xDiag(end) );
yDiag(end) = ( 0.4*yDiag(end-1)+0.6*yDiag(end) );

% Combine teeth parts into one single vector
xTooth = [xH xDiag];
yTooth = [yH yDiag];

%
% TESTING GEOMETRY
%plot(xH,yH,'r.'); hold on;
%plot(xDiag,yDiag,'r.'); hold on;
%plot(xTooth,yTooth,'k'); hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: rotate the geometry by ang
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xDiagN,yDiagN] = rotate_geometry(ang,xDiag,yDiag)

%
% Rotate Diagonal Piece
for i=1:length(xDiag)
   xDiagN(i) = xDiag(i)*cos(ang) - yDiag(i)*sin(ang);
   yDiagN(i) = xDiag(i)*sin(ang) + yDiag(i)*cos(ang);
end

