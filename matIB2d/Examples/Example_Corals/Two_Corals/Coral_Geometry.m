function Coral_Geometry()

L = 0.5;
Nx = 128;
Ny = 128;
ratio = 1024 / (2*Nx);

struct_name = 'coral'; % Name for .vertex, .spring, etc files.

[x1,y1,N] = please_Give_Single_Polyp_Geometry(ratio);
x1=x1-0.1; x2=x1+0.2;
y2=y1;
x = [x1;x2];
y = [y1;y2];

fprintf('\n\nNumber of Pts. in ONE tentacle: %d\n\n',N);

% Prints .vertex file!
print_Lagrangian_Vertices(x,y,struct_name);

% Prints .target file!
k_Target = 1e8;
print_Lagrangian_Target_Pts(x,k_Target,struct_name);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints VERTEX points to a file called corals.vertex
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
% FUNCTION: prints Vertex points to a file called corals.vertex
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in the vertices for one coral polyp
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,y,N] = please_Give_Single_Polyp_Geometry(ratio)

xoffset = 0.3; % To put geometry into QUADRANT-1
yoffset = 0.5; % To put geometry into QUADRANT-1

% Note: (1) ds for stored data is 0.6/(2*1024)
%       (2) 'ratio' is comparing 1024:desired resolution

struct_name1 = 'coral2d_left_1024';
[~,x1,y1] = read_Vertex_Points(struct_name1);

struct_name2 = 'coral2d_right_1024';
[~,x2,y2] = read_Vertex_Points(struct_name2);

% Put Geometry Together for One Polyp
xAux = [x1+xoffset; x2(end:-1:1)+xoffset];
yAux = [y1+yoffset; y2(end:-1:1)+yoffset];

x=xAux(1:ratio:end);
y=yAux(1:ratio:end);
N = length(x)/2;

plot(xAux,yAux,'*'); hold on;
plot(x,y,'ro'); hold on;

%NOTE: (1) N here is the # of pts. on one arm (NOT ENTIRE POLYP)!!
%      (2) The geometry here is listed for 1024x1024 meshes -> will need
%          to take every 8 or 16 pts to render geometry usable for MATLAB
%          code


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in the # of vertex pts and all the vertex pts from the
%           .vertex file.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [N,xLag,yLag] = read_Vertex_Points(struct_name)

filename = [struct_name '.vertex'];  %Name of file to read in
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