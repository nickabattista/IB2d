function Coral_Geometry()

Lx = 2;
Ly = 5;
Nx = 256;
Ny = 640;

dt=2.5e-4;
tf=18.9;
t=dt:dt:tf;
dx=Lx/Nx;
dy=Ly/Ny;

struct_name = 'coral'; % Name for .vertex, .spring, etc files.

[x,y,N] = please_Give_Single_Polyp_Geometry(dx,dt,t);


fprintf('\n\nNumber of Pts. in ONE tentacle: %d\n\n',N);

% Prints .vertex file!
print_Lagrangian_Vertices(x,y,struct_name);

% Prints .target file!
k_Target = 8e5;
print_Lagrangian_Target_Pts(x,k_Target,struct_name);

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

xoffset = 1; % To put geometry into QUADRANT-1
yoffset = 0.5; % To put geometry into QUADRANT-1

% Note: (1) ds for stored data is 0.6/(2*1024)
%       (2) 'ratio' is comparing 1024:desired resolution

ds=h/2;

load('total_coeffs.mat')
load('coral_coeff_30.mat')
t=t/total_pulse;


%arclength
s=(0:ds:total_meanL)/total_meanL;

% number of points on one arm
N=length(s)

c1vals=[];
c2vals=[];
tavgex=[t t(end)+dt];
tavg=(tavgex(2:end)+tavgex(1:end-1))/2;
tavg=[0 tavg];
size(cs1(:,1))
for indx1=1:4
c1vals(:,indx1)=ppval(cs1(:,indx1),tavg);
c2vals(:,indx1)=ppval(cs2(:,indx1),tavg);
end

C1=c1vals(1,:);
C2=c2vals(1,:);

     XbL_1=(C1(1)*s.^3+C1(2)*s.^2+C1(3)*s+C1(4));
     YbL_1=(C2(1)*s.^3+C2(2)*s.^2+C2(3)*s+C2(4));

     L_ten=sum(sqrt((XbL_1(2:end)-XbL_1(1:end-1)).^2 +(YbL_1(2:end)-YbL_1(1:end-1)).^2 ));

     XbL=(C1(1)*s.^3+C1(2)*s.^2+C1(3)*s+C1(4))*total_meanL/L_ten+total_offset;
     YbL=(C2(1)*s.^3+C2(2)*s.^2+C2(3)*s+C2(4))*total_meanL/L_ten;

    XbR=-XbL;
    YbR=YbL;


XbStem=(XbR(1)+ds*3/2):ds:(XbL(1)-ds);
YbStem=(YbL(1))+0*XbStem;

Nstem=length(XbStem);
save('cval.mat','c1vals','c2vals','Nstem','XbStem','YbStem')

x1=[flip(XbR) XbStem XbL];
y1=[flip(YbR) YbStem YbL];

% Put Geometry Together for One Polyp
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
% FUNCTION: gives initial concentration gradient inside channel
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C,X,Y] = give_Me_Initial_Concentration(Lx,Ly,Nx,Ny,dx,dy,xLag,yLag)
%X_lag(:,1)=xLag(1:ceil(length(xLag)/2));
%X_lag(:,2)=xLag(ceil(length(xLag)/2)+1:end);

%Y_lag(:,1)=yLag(1:ceil(length(xLag)/2));
%Y_lag(:,2)=yLag(ceil(length(xLag)/2)+1:end);


x = 0:dx:Lx-dx;
y = 0:dy:Ly-dy;

[X,Y]=meshgrid(x,y);

C=0*X;
%eps=2*dx;
%f=X_lag(:,1)*0+1/20;

%for ten=1:2
%XL=X_lag(:,ten);
%YL=Y_lag(:,ten);

%for i=1:Nx
%    for j=1:Nx
%        %find x and y
%
%        xi=X(i,j);
%        yj=Y(i,j);
%
%        % regularized delta function
%        delta=1/((eps*sqrt(2*pi)).^2)*exp(-(xi-XL).^2/(2*eps^2)).*exp(-(yj-YL).^2/(2*eps^2));
%        w=f.*delta;
%        
%        ds_vec=cumsum(sqrt((XL(2:end)-XL(1:end-1)).^2+(YL(2:end)-YL(1:end-1)).^2));
%        
%        ds_vec=[0; ds_vec];
%        
%        C_test(i,j,ten)=trapz(ds_vec,w);
%
%        %f(i,j)=trapz(theta,F.*w);
%        %f2(i,j)=(1/2*w(1).*F2(1)+sum(w(2:end-1).*F2(2:end-1))+1/2*w(end).*F2(end))*dtheta;
%
%    end
%end
%end 
%C=C_test(:,:,1)+C_test(:,:,2);
