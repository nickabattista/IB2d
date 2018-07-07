%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in the # of vertex pts and all the vertex pts from the
%           .vertex file and plots them.
%
% NOTE: pass string of geometry name in as input, e.g., 'rubberband' to
%       read in and plot rubberband.vertex points.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Leafitize_Me()

% Note: Leaf Geometry Coming in at Resolution 256x256
N = 1024;
L = 2;
ds = (L)/(2*N);
struct_name = 'leaf_vein';

%
% Lagrangian Pts at 128x256 Resolution on [0,1]x[0,2]
%
[xLag,yLag] = read_Vertex_Points_and_Plot_Them('leaf_flow');


%
% Make higher resolution geometry (always do factor of 2...then nest, if
% necessary)
%
[xLag,yLag] = make_Higher_Resolution(xLag,yLag,2*ds);
[xLag,yLag] = make_Higher_Resolution(xLag,yLag,ds);

%
% Add tiny channel at bottom for inflow
%
[xLag,yLag] = tack_On_Vertical_Channel_Lines(xLag,yLag,N,L,ds);

%
% Prints .vertex file!
%
print_Lagrangian_Vertices(xLag,yLag,struct_name);

%
% Prints .target file!
%
k_Target = 1e5;
print_Lagrangian_Target_Pts(xLag,k_Target,struct_name);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: adds Vertical Channel Lines onto the leaf vein geometry 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag] = make_Higher_Resolution(xLag,yLag,ds)
   
    % ALWAYS RESOLVES BY A FACTOR OF 2!

    xLag_N = zeros(1,2*length(xLag)-1);
    yLag_N = xLag_N;
    
    %
    % Put pts in the middle between successive points
    %
    for i=1:length(xLag)-1

        xCur = xLag(i);
        yCur = yLag(i);

        xNext = xLag(i+1);
        yNext = yLag(i+1);

        xMid = ( xCur + xNext ) / 2;
        yMid = ( yCur + yNext ) / 2;

        dist = sqrt( (xMid-xCur)^2 + (yMid-yCur)^2 );

        if dist < 2*ds
            xLag_N(2*i-1) = xLag(i);
            yLag_N(2*i-1) = yLag(i);

            xLag_N(2*i) = xMid; 
            yLag_N(2*i) = yMid;
        else

            xLag_N(2*i-1) = xLag(i);
            yLag_N(2*i-1) = yLag(i);

            xLag_N(2*i) = 5000; % dummy value
            yLag_N(2*i) = 5000; % dummy value

        end  

    end

    % Remove all entries with 5000 in it
    xLag = xLag(xLag~=5000);
    yLag = yLag(yLag~=5000);
    
    figure(2)
    plot(xLag,yLag,'.'); hold on
    
    % check
    %Inds = find(xLag==5000)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: adds Vertical Channel Lines onto the leaf vein geometry 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag] = tack_On_Vertical_Channel_Lines(xLag,yLag,N,L,ds)

%yLow = yLag(1);
%for i=2:length(xLag)
%    if yLag(i)< yLow
%        yLow = yLag(i);
%    end
%end

yLow = 0.219;
n = 1;
for i=1:length(xLag)
    if yLag(i) > yLow
        xLagN(n) = xLag(i);
        yLagN(n) = yLag(i);
        n=n+1;
    end
end


yAdd=0:ds:L/10;
xAdd = ones(1,length(yAdd));

% Left / Right Sides for inflow channel
xLeft = 0.247375;
xRight = 0.35955;

yAdds = [yAdd yAdd];
yAdds = ( yLow-(L/10)-ds/2 )+ yAdds;
xAdds = [xLeft*xAdd xRight*xAdd];

xLag = [xLagN xAdds];
yLag = [yLagN yAdds];

yLag = yLag - L/20;

% Center in Window [0,1]x[0,2.5]
yLag = yLag + L/8;

figure(1)
%plot(xLagN,yLagN,'g*'); hold on;
%plot(xAdds,yAdds,'k*'); hold on;
plot(xLag,yLag,'r*'); hold on;
%plot(xLag,yLag,'bo'); hold on;
axis([0 1 0 2.5]);

fprintf('\n\n        FOR SETTING UP INFLOW:\n\n');
fprintf('   (x,y)-LEFT: (%d,0.15)\n',xLeft);
fprintf('   (x,y)-RIGHT: (%d,0.15)\n\n\n',xRight);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: reads in vertex points from "leaf_flow.vertex" and passes them
% back to modify them. NOTE: original resolution at 256x256
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag] = read_Vertex_Points_and_Plot_Them(struct_name)

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