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

function Jellyfish_Geometry()

% FLUID GRID PARAMETERS %
L = 8;                              % Length of computational domain (m)
N = 512;                            % # of Cartesian grid meshwidths
dx = L/N;                           % Cartesian mesh width (m)

% Construct Geometry
[xLag,yLag,ds] = give_Me_Immsersed_Boundary_Geometry(N,L);

for i =1:length(xLag)
    plot(xLag(i),yLag(i),'*'); hold on;
    pause(0.1);
end

% Translate Geometry
xLag = xLag + L/8;
yLag = yLag + L/4;
plot(xLag,yLag,'b*'); hold on;
plot(xLag(1),yLag(1),'m*'); hold on;
plot(xLag(127),yLag(127),'y*'); hold on;
plot(xLag(end),yLag(end),'g*'); hold on;

  
% NAMING CONVENTION FOR SIMULATION 
struct_name = 'jelly';      % structure name


%
% PRINT INPUT FILES (.vertex, .spring, .beam, etc) %
%

% print vertices
print_Lagrangian_Vertices(xLag,yLag,struct_name);

% print springs
k_Spring = 1.2750000000000000e+07;   % spring constant (Newton)
print_Lagrangian_Springs(xLag,k_Spring,ds,struct_name);

% print beams
k_Beam = 1.0363359375000002e+15;   % beam stiffness constant (Newton m^2)
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
% FUNCTION: prints SPRING points to a file called rubberband.spring
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Springs(xLag,k_Spring,ds_Rest,struct_name)

    N = length(xLag);  % N IS ODD

    spring_fid = fopen([struct_name '.spring'], 'w');

    fprintf(spring_fid, '%d\n', N-1 );

    %spring_force = kappa_spring*ds/(ds^2);

    % SPRINGS BETWEEN VERTICES ON RHS
    for s = 1:(N-1)/2
            %if s < (N-1)/2+1        
                fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, s+1, k_Spring, ds_Rest);  
            %elseif s == ceil(N/2)
            %    fprintf(spring_fid, '%d %d %1.16e %1.16e\n', 1, ceil(N/2)+1, k_Spring, ds_Rest);  
            %end
    end
    
    % SPRINGS BETWEEN VERTICES ON LHS
    for s=(N-1)/2+2:N
            
        if s==(N-1)/2+2
            s1 = 1;%ceil(N/2)+s;
            s2 = s;%ceil(N/2)+s+1;
            fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s1, s2, k_Spring, ds_Rest);  
        else
            s1 = s-1;%ceil(N/2)+s;
            s2 = s;%ceil(N/2)+s+1;
            fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s1, s2, k_Spring, ds_Rest);      
        end
    end
    fclose(spring_fid); 
    
 
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
    for s = 1:(N-1)/2
        if s==1
            Cx = xLag((N-1)/2+2) - 2*xLag(1) + xLag(2);
            Cy = yLag((N-1)/2+2) - 2*yLag(1) + yLag(2);
            fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', (N-1)/2+2,1,2, k_Beam, Cx,Cy);
        elseif s <= (N-1)/2  
            Cx = xLag(s-1) - 2*xLag(s) + xLag(s+1);
            Cy = yLag(s-1) - 2*yLag(s) + yLag(s+1);
            fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', s-1,s,s+1, k_Beam, Cx, Cy); 
        end
    end
    
    % BEAMS BETWEEN VERTICES ON LHS
    for s=(N-1)/2+2:N-1
         if s==(N-1)/2+2
            s1 = 1;
            s2 = s;
            s3 = s+1;
            Cx = xLag(s3) - 2*xLag(s2) + xLag(s1);
            Cy = yLag(s3) - 2*yLag(s2) + yLag(s1);
            fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', s3, s2, s1, k_Beam, Cx, Cy);  
         else
            s1 = s-1;
            s2 = s;
            s3 = s+1;
            Cx = xLag(s3) - 2*xLag(s2) + xLag(s1);
            Cy = yLag(s3) - 2*yLag(s2) + yLag(s1);
            fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', s3, s2, s1, k_Beam, Cx, Cy);  
         end
    end
    fclose(beam_fid); 
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates the Lagrangian structure geometry
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag,ds] = give_Me_Immsersed_Boundary_Geometry(N,L)


% JELLYFISH GEOMETRY PARAMETERS %

bell_length = 2;                    % bell length (m)
npts_bell = ceil(2*(bell_length/L)*N); % number of points along the length of the entire bell
npts_circ = 1;                    %number of points along the circumference (=1 for 2D)
npts = npts_bell*npts_circ;	      % total number points 
ds = bell_length/(npts_bell-1);   % mesh spacing along the length of bell (m)
Zs = 0;                          %distance to top of bell (m)
xb = zeros(npts);                 %holds resting positions to calculate curvatures
zb=zeros(npts);                   %holds resting positions to calculate curvatures

% Values to describe bell shape given in Alben, Peng and Miller %
betao = 0.5;
betam = 0.3;
to = 0.5;
t=0;
gamma = 1;

% Starting values
z = Zs;
r = 0;
x = 0;

% center line of jelly
xLag(1) = x;
yLag(1) = z;
zl=z;
rl=r;

%right side of bell
for s = 1:(ceil(npts_bell/2)-1)
    beta = betao+(betam-betao)*(t/to)^gamma;
    theta = -1.55*(1-exp(-(s)*ds/beta));
    z = zl + ds*sin(theta);
    r = rl + ds*cos(theta);
    x = r;
    zl=z;
    rl=r;
    xLag(s+1)=x;
    yLag(s+1)=z;
   %fprintf(vertex_fid, '%1.16e %1.16e\n', x, z);
end


% reinitialize values for centerline
z = Zs;
r = 0;
zl=z;
rl=r;


%left side of bell
for s = (ceil(npts_bell/2)):npts_bell-2
    s2=s-(ceil(npts_bell/2)-1);
    beta = betao+(betam-betao)*(t/to)^gamma;
    theta = -1.55*(1-exp(-(s2)*ds/beta));
    z = zl + ds*sin(theta);
    r = rl + ds*cos(theta);
    x = -1*r;
    zl=z;
    rl=r;
    xLag(s+1)=x;
    yLag(s+1)=z;
   %fprintf(vertex_fid, '%1.16e %1.16e\n', x, z);
end


% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Write out the beam information
%  beam_fid = fopen([mesh_name num2str(N) '.beam'], 'w');
%  
%  fprintf(beam_fid, '%d\n', npts-3);
%  
%  %right side of bell
%  for q = 0:npts_circ-1
%  for s = 0:(ceil(npts_bell/2)-3)
%     C1 = xb(s+1)+xb(s+3)-2*xb(s+2); 
%     C2 = zb(s+1)+zb(s+3)-2*zb(s+2); 
%     fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', q*npts_circ+s, q*npts_circ+s+1, q*npts_circ+s+2, kappa_beam*ds/(ds^4), C1, C2);
%  end
%  
%  %top of bell
%     s=ceil(npts_bell/2);
%     C1 = xb(s+1)+xb(2)-2*xb(1); 
%     C2 = zb(s+1)+zb(2)-2*zb(1); 
%     fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', q*npts_circ+s, q*npts_circ+0, q*npts_circ+1, kappa_beam*ds/(ds^4), C1, C2);
%     s=ceil(npts_bell/2)+1;
%     C1 = xb(s+1)+xb(1)-2*xb(s); 
%     C2 = zb(s+1)+zb(1)-2*zb(s); 
%     fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', q*npts_circ+s, q*npts_circ+s-1, q*npts_circ+0, kappa_beam*ds/(ds^4), C1, C2);
%  
%  %left side of bell   
%  for s = ceil((npts_bell/2)+2):npts_bell-2
%     C1 = xb(s+1)+xb(s-1)-2*xb(s); 
%     C2 = zb(s+1)+zb(s-1)-2*zb(s); 
%     fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', q*npts_circ+s, q*npts_circ+s-1, q*npts_circ+s-2, kappa_beam*ds/(ds^4), C1, C2);
%  end
%  end
% % 
% % 
%  fclose(beam_fid);
% 
