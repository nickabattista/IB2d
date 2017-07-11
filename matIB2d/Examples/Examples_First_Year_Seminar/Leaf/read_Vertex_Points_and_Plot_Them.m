%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in the # of vertex pts and all the vertex pts from the
%           .vertex file and plots them.
%
% NOTE: pass string of geometry name in as input, e.g., 'rubberband' to
%       read in and plot rubberband.vertex points.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  read_Vertex_Points_and_Plot_Them(struct_name)

%struct_name = (name of geometry), e.g., rubberband for rubberband.vertex

filename = [struct_name '.vertex_OLD'];  %Name of file to read in
fileID = fopen(filename);

% Read in the file, use 'CollectOutput' to gather all similar data together
% and 'CommentStyle' to to end and be able to skip lines in file.
C = textscan(fileID,'%f %f','CollectOutput',1);

fclose(fileID);     %Close the data file.

vertices = C{1};    %Stores all read in data in vertices (N+1,2) array

N = vertices(1,1);  % # of Lagrangian Pts
xLag = zeros(N,1);  % Initialize storage for Lagrangian Pts.
yLag = xLag;        % Initialize storage for Lagrangian Pts.

for i=1:N
   xLag(i,1) = vertices(i+1,1); %Stores x-values of Lagrangian Mesh
   yLag(i,1) = vertices(i+1,2); %Stores y-values of Lagrangian Mesh
end

% Center the geometry in domain and push it up a bit
xLag = xLag + 0.0125;
yLag = yLag + 0.035;

% Rotate the geometry
[xLag,yLag] = please_Rotate_Geometry(xLag,yLag);

xLag = xLag + 0.05;
yLag = yLag + 0.03; % just to push it up a bit further in domain :)

plot(xLag,yLag,'r*'); hold on;
plot(xLag,yLag,'bo'); hold on;
axis([0 0.2 0 0.3]); hold on;

print_Lagrangian_Vertices(xLag,yLag,struct_name)

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
% FUNCTION: rotates the geometry by an angle, ang.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [xLag,yLag] = please_Rotate_Geometry(xLag,yLag)
        
% ROTATE GEOMETRY 180 DEGREES: approximate center of geometry is (xC,yC) = (0.05, 0.075)
ang = pi;   % angle for rotation
xC = 0.05;  % x-Center of rotation
yC = 0.075; % y-Center of rotation

% PERFORM ROTATION
xLagN = (xLag-xC)*cos(ang) - (yLag-yC)*sin(ang);
yLagN = (xLag-xC)*sin(ang) + (yLag-yC)*cos(ang);

% PUT BACK IN CENTER OF DOMAIN
xLag = xLagN+xC;
yLag = yLagN+yC;