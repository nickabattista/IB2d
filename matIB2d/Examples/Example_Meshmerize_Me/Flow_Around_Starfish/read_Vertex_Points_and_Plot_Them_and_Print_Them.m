%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: reads in vertex points from "star_fish.vertex" and passes them
%           back to modify them. 
%
%           NOTE: original resolution set at 256x256 for star_flow.vertex
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function read_Vertex_Points_and_Plot_Them_and_Print_Them(struct_name)

% 
% Read in Lagrangian Pts from Given Vertex File
%
[xLag,yLag] = retrieve_Lagrangian_Pts_From_Vertex(struct_name);

%
% Lengths of (Desired) Computational Domain / Grid Resolution 
%
Lx = 1.0; Ly = 1.0;
Nx = 1024; Ny = 1024;

%
% Proxy for Height of Object (scale factor)
% 
hx = 0.25;  % Object is supposed to take up "approximately 1/4th region of original domain (L=1)"

%
% Desired Lagrangian Spacing, ds
%
ds_Desired = Lx/(2*Nx);

%
% Scales the Points (by a scale factor) that is determined by taking the
% original "equally" spaced pts and scaling them appropriately to give the
% desired Lagrangian spacing set above.
%
ds_Orig = mean (sqrt( (xLag(1:end-1) - xLag(2:end)).^2 + (yLag(1:end-1)-yLag(2:end)).^2  ) );
scale = ds_Desired/ds_Orig;
xLagN = scale*xLag; yLagN = scale*yLag;
ds_Check_Vec =  mean (sqrt( (xLagN(1:end-1) - xLagN(2:end)).^2 + (yLagN(1:end-1)-yLagN(2:end)).^2  ) );
ds_Check = mean( ds_Check_Vec );
ds_err = abs(ds_Check-ds_Desired);

fprintf('\n\nOriginal Lagrangian Spacing: %d\n',ds_Orig);
fprintf('Desired Lagrangian Spacing: %d\n',ds_Desired);
fprintf('Scale Factor: %d\n',scale);
fprintf('New Lagrangian Spacing, ds (check): %d\n',ds_Check);
fprintf('Lagrangian Spacing Avg-Error: %d\n',ds_err);
fprintf('Lagrangian Spacing Max-Error: %d\n\n\n',max(ds_Check_Vec-ds_Desired));


%
% Finds Center Pt. of Starfish for Easier Translation of Geometry
%
xC = mean(xLagN);
yC = mean(yLagN);

%
% Translates Points to Middle of Domain
%
xLagN = xLagN + (Lx/2-xC) + Lx/2; 
yLagN = yLagN + (Ly/2-yC) - Ly/4;



%
% Plot to Test Geometry
%
figure(1)
plot(xLag,yLag,'b.-'); hold on;
plot(xLagN,yLagN,'r.-'); hold on;
legend('Original Size','New Size');
xlabel('x');
ylabel('y');
%axis([0 Lx 0 Ly])
axis([0 4 0 0.5]);
%axis square;

%
% PRINTS NEW .VERTEX and .TARGET FILES
%
struct_name = 'starfish_mod';
k_Target = 1e7;
print_Lagrangian_Vertices(xLagN,yLagN,struct_name);
print_Lagrangian_Target_Pts(xLag,k_Target,struct_name);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: get all Lagrangian Pts from Original .vertex file {(xLag_j,yLag_j)}_j=1^N 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag] = retrieve_Lagrangian_Pts_From_Vertex(struct_name)

%struct_name = (name of geometry), e.g., rubberband for rubberband.vertex

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

% Store vertices in appropriate xLag, yLag arrays
for i=1:N
   xLag(i,1) = vertices(i+1,1); %Stores x-values of Lagrangian Mesh
   yLag(i,1) = vertices(i+1,2); %Stores y-values of Lagrangian Mesh
end


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
% FUNCTION: prints Vertex points to a file called <struct_name>.target
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
