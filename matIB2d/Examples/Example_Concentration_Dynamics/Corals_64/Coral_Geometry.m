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
% IB2d Date Created: May 27th, 2015
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates the input files (.vertex, .target, etc.) for the 
%           Coral polyp example
%
%   Authors: Dr. Matea Santiago, Dr. Shilpa Khatri
%   Created: March 2021
%   Modified: July 2021 (by Dr. Nick Battista)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Coral_Geometry()

Lx = 2;    % Length of horizontal domain
Ly = 5;    % Length of vertical domain 
Nx = 64;  % Resolution in horizontal direction
Ny = 160;  % Resolution in vertical direction
dx=Lx/Nx;  % Grid spacing in horizontal direction
dy=Ly/Ny;  % Grid spacing in vertical direction

dt=2.5e-4;  % Time-step
tf=18.9;    % Final time in simulation
t=dt:dt:tf; % Vector of times

%----------------------------------
% Structure name for input files
%----------------------------------
struct_name = 'coral'; % Name for .vertex, .spring, etc files.


%---------------------------------------------------------
% Give x and y positions for 1st state and # of points
%   on ONE tentacle
%---------------------------------------------------------
[x,y,N] = please_Give_Single_Polyp_Geometry(dx,dt,t);


fprintf('\nNumber of Pts. on ONE tentacle: %d\n',N);
fprintf('\nNumber of Pts. in full CORAL POLYP: %d\n\n',length(x));
fprintf('---> Make sure the UPDATE_TARGET_POINT_POSITIONS.m\n');
fprintf('     has the correct Nx and Lx\n');
fprintf('    (the same as in this file and input2d)\n\n');


% Prints .vertex file!
print_Lagrangian_Vertices(x,y,struct_name);

% Prints .target file!
k_Target = 8e5;
print_Lagrangian_Target_Pts(x,k_Target,struct_name);

% Prints Geometry Connections!
print_Geometry_Connections(x,struct_name);

% Call function for initial concentration
[Concentration,X,Y]= give_Me_Initial_Concentration(Lx,Ly,Nx,Ny,dx,dy,x,y);

% Prints .concentration file!
%kDiffusion = 1e-2;
print_Concentration_Info(Nx,Ny,Concentration,struct_name);





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

function [x,y,N] = please_Give_Single_Polyp_Geometry(h,dt,t)

xoffset = 1;   % To put geometry into QUADRANT-1
yoffset = 0.5; % To put geometry into QUADRANT-1
ds=h/2;        % Lagrangian point spacing

%-----------------------------------------------
% Note: (1) ds for stored data is 0.6/(2*1024)
%       (2) 'ratio' is comparing 1024:desired resolution
%-----------------------------------------------

%-----------------------------------------------
% Load prescribed position data for tentacles
%-----------------------------------------------
load('total_coeffs.mat')
load('coral_coeff_30.mat')
t=t/total_pulse;

%-----------------------------------------------
% Arclength
%-----------------------------------------------
s=(0:ds:total_meanL)/total_meanL;

%-----------------------------------------------
% Number of points on one tentacle
%-----------------------------------------------
N=length(s);


%-----------------------------------------------
% Set up interpolate state times
%-----------------------------------------------
c1vals=[];
c2vals=[];
tavgex=[t t(end)+dt];
tavg=(tavgex(2:end)+tavgex(1:end-1))/2;
tavg=[0 tavg];

%----------------------------------------------------------------
% Get interpolation polynomial coefficients to interpolate btwn
%      tentacle geometries
%-----------------------------------------------------------------
for indx1=1:4
    c1vals(:,indx1)=ppval(cs1(:,indx1),tavg);
    c2vals(:,indx1)=ppval(cs2(:,indx1),tavg);
end
C1=c1vals(1,:);
C2=c2vals(1,:);

%-----------------------------------------------
% Get LEFT tentacle interpolation states
%-----------------------------------------------
XbL_1=(C1(1)*s.^3+C1(2)*s.^2+C1(3)*s+C1(4));
YbL_1=(C2(1)*s.^3+C2(2)*s.^2+C2(3)*s+C2(4));

L_ten=sum(sqrt((XbL_1(2:end)-XbL_1(1:end-1)).^2 +(YbL_1(2:end)-YbL_1(1:end-1)).^2 ));

XbL=(C1(1)*s.^3+C1(2)*s.^2+C1(3)*s+C1(4))*total_meanL/L_ten+total_offset;
YbL=(C2(1)*s.^3+C2(2)*s.^2+C2(3)*s+C2(4))*total_meanL/L_ten;

%-------------------------------------------------------
% Set up symmetric RIGHT tentacle interpolation states
%-------------------------------------------------------
XbR=-XbL;
YbR=YbL;

%-------------------------------------------------------
% Setup polyp STEM geometry
%-------------------------------------------------------
XbStem=(XbR(1)+ds*3/2):ds:(XbL(1)-ds);
YbStem=(YbL(1))+0*XbStem;
Nstem=length(XbStem);

%-------------------------------------------------------
% Save positional state data to .mat files
%-------------------------------------------------------
save('cval.mat','c1vals','c2vals','Nstem','XbStem','YbStem')

%-------------------------------------------------------
% Setup FIRST state for coral polyp to begin simulation
%-------------------------------------------------------
x1=[flip(XbR) XbStem XbL];
y1=[flip(YbR) YbStem YbL];

%-------------------------------------------------------
% Put Geometry Together for One Polyp and translate 
%       to desired place in domain appropriately
%-------------------------------------------------------
x = x1+xoffset;
y = y1+yoffset;
plot(x,y,'.')
axis([0 2 0 2])
drawnow

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints CONCENTRATION INFO to file called
%           'struct_name'.concentration
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function print_Concentration_Info(Nx,Ny,C,struct_name)

    con_fid = fopen([struct_name '.concentration'], 'w');

%    fprintf(con_fid, '%d\n', kDiffusion );

    for i=1:Ny
        for j=1:Nx
            fprintf(con_fid, '%1.16e ', C(i,j) );
        end
        fprintf(con_fid,'\n');
    end    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: print geometrical connections between adjacent 
%           Lagrangian points
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function print_Geometry_Connections(xLag,struct_name)


    N = length(xLag); % Total # of Lag. Pts

    geo_fid = fopen([struct_name '.geo_connect'], 'w');

    % Loops over all Lagrangian Pts.
    for i = 1:N
        
        if i<N
            
            s1 = i;   
            s2 = i+1;
            fprintf(geo_fid, '%d %d\n', s1, s2);
            fprintf(geo_fid, '%d %d\n', s2, s1);
            
        end
        
    end


    fclose(geo_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives initial concentration in domain
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C,X,Y] = give_Me_Initial_Concentration(Lx,Ly,Nx,Ny,dx,dy,xLag,yLag)

x = 0:dx:Lx-dx;
y = 0:dy:Ly-dy;

[X,Y]=meshgrid(x,y);

C=zeros(Ny,Nx);
