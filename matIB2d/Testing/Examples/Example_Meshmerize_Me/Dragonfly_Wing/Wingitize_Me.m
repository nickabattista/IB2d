%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in the # of vertex pts and all the vertex pts from the
%           .vertex file and plots them.
%
% NOTE: pass string of geometry name in as input, e.g., 'rubberband' to
%       read in and plot rubberband.vertex points.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Wingitize_Me()

% Note: Wing Geometry Coming in at Resolution 1024x410 on [0,5]x[0,2]
N = 1024;
Lx = 0.5;
Ly = 0.5;
ds = (Lx)/(2*N);
struct_name = 'wing_vein';

%
% Lagrangian Pts at 128x256 Resolution on [0,1]x[0,2]
%
[xLag,yLag] = read_Vertex_Points_and_Plot_Them('wing_flow');

%
% Get rid of geometry outside of range!
%
xLeft = 0.0045; yLow = 0.01;
[xLag,yLag] = get_Rid_Of_Erroneous_Points(yLow,xLeft,xLag,yLag);

xLag = xLag+0.02;
yLag = yLag+0.02;


figure(3)
plot(xLag,yLag,'.'); hold on;
axis([0 .5 0 .25])
%pause();

for i=1:length(xLag)-1
   dsVec(i) =  sqrt( (xLag(i)-xLag(i+1))^2 + (yLag(i)-yLag(i+1))^2   );
end

% original ds from MeshmerizeMe
ds = median(dsVec);


%
% Make higher resolution geometry (always do factor of 2...then nest, if
% necessary)
%
[xLag,yLag] = make_Higher_Resolution(xLag,yLag,ds);
[xLag,yLag] = make_Higher_Resolution(xLag,yLag,ds);
%[xLag,yLag] = make_Higher_Resolution(xLag,yLag,ds);

for i=1:length(xLag)-1
   dsVec2(i) =  sqrt( (xLag(i)-xLag(i+1))^2 + (yLag(i)-yLag(i+1))^2   );
end

% Modified ds to new geometry
ds = median(dsVec2);



%
% Add tiny channel at bottom for inflow
%
[xLag,yLag] = tack_On_Vertical_Channel_Lines(xLag,yLag,N,Lx,ds);

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
% FUNCTION: gets rid of erroneous points for insect wing veins 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag] = get_Rid_Of_Erroneous_Points(yLow,xLeft,xLag,yLag)

%
% GET RID OF POINTS W/ y<0:
%
ct = 1;
for i=1:length(xLag)
    if yLag(i) > yLow 
       xLagN(ct) = xLag(i);
       yLagN(ct) = yLag(i);
       ct = ct+1;
    end
end

% Redefine Lag. Pts
xLag = xLagN; yLag = yLagN;
%plot(xLagN,yLagN,'r.'); hold on;

%
% GET RID OF POINTS TO LEFT OF X=0.004 (for beginning of channel off end)
%
ct = 1;
for i=1:length(xLag)
    if xLag(i) > xLeft 
       xLagN(ct) = xLag(i);
       yLagN(ct) = yLag(i);
       ct = ct+1;
    end
end

xLag = xLagN;
yLag = yLagN;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: adds Vertical Channel Lines onto the leaf vein geometry 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag] = make_Higher_Resolution(xLag,yLag,ds)
   
    % ALWAYS RESOLVES BY A FACTOR OF 2!

    xLag_N = zeros(1,2*length(xLag)-2);
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

            xLag_N(2*i) = 5000.0; % dummy value
            yLag_N(2*i) = 5000.0; % dummy value

        end  

    end

    % Remove all entries with 5000 in it
    xLag = xLag_N(xLag_N~=5000.0);
    yLag = yLag_N(yLag_N~=5000.0);
    
    %figure(2)
    %plot(xLag,yLag,'r.'); hold on
    
    ct=0;
    for i=1:length(xLag)
        if xLag(i)==0
            ct=ct+1;
        end
    end
    
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

figure(4)
plot(xLag,yLag,'ro'); hold on;

xRight = 0.02465;
yBotVal = 0.19210; 
yTopVal = 0.19827;

xTop = xRight:-ds:0.01;
xBot = xTop;

yBot = yBotVal*ones(1,length(xTop));
yTop = yTopVal*ones(1,length(xTop));

xLag = [xLag xBot xTop];
yLag = [yLag yBot yTop];

plot(xLag,yLag,'g.'); hold on;


%figure(1)
%plot(xLagN,yLagN,'g*'); hold on;
%plot(xAdds,yAdds,'k*'); hold on;
%plot(xLag,yLag,'r*'); hold on;
%plot(xLag,yLag,'bo'); hold on;
%axis([0 1 0 2.5]);

fprintf('\n\n        FOR SETTING UP INFLOW:\n\n');
fprintf('   (x,y)-TOP: (%d,0.025)\n',yTopVal);
fprintf('   (x,y)-BOTTOM: (%d,0.025)\n\n\n',yBotVal);



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