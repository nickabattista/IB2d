%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  nick.battista@unc.edu
% Date Created: September 9th, 2016
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
% If you would like us %to add a specific muscle model, please let Nick (nick.battista@unc.edu) know.
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates the geometry and prints associated input files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Make_Tub_Geometry_and_Initial_Concentration()

%
% Grid Parameters (MAKE SURE MATCHES IN input2d !!!)
%
Nx =  128;       % # of Eulerian Grid Pts. in x-Direction (MUST BE EVEN!!!)
Lx = 1.0;        % Length of Eulerian Grid in x-Direction
dx = Lx/Nx;      % Grid spatial resolution

%
% Immersed Structure Geometric / Dynamic Parameters %
%
ds = 0.25*dx;              % Lagrangian Pt. Spacing (2x resolution of Eulerian grid)
struct_name = 'boussinesq'; % Name for .vertex, .spring, etc files. (must match what's in 'input2d')


% Call function to construct geometry
[xLag,yLag,xSpider,ySpider,Concentration,Nspider] = give_Me_Immsersed_Boundary_Geometry_and_Concentration(ds,Nx,Lx,dx);



% Plot Geometry to test
plot(xLag,yLag,'r-'); hold on;
plot(xLag,yLag,'*'); hold on;
xlabel('x'); ylabel('y');
axis square;



% Prints .vertex file!
print_Lagrangian_Vertices([xLag xSpider],[yLag ySpider],struct_name);


% Prints .spring file!
k_Spring = 2.5e6;                   % Spring stiffness (does not need to be equal for all springs)
ds_Rest = ds;                       % Spring resting length (does not need to be equal for all springs)
print_Lagrangian_Springs([xLag xSpider],yLag,k_Spring,ds_Rest,struct_name,length(xLag),Nspider);


% Prints .beam file!
k_Beam = 2e5;                      % Beam Stiffness (does not need to be equal for all beams)
%C = compute_Curvatures(xLag,yLag)  % Computes curvature of initial configuration
C = 0;
print_Lagrangian_Beams([xLag xSpider],yLag,k_Beam,C,struct_name,length(xLag),Nspider);

% Prints .mass file!
k_Mass = 1e5; Mass = 2e-4;
print_Lagrangian_Mass_Pts([xLag xSpider],k_Mass,Mass,struct_name,length(xLag),Nspider)


% Prints .target file!
k_Target = 1e7;
print_Lagrangian_Target_Pts(xLag,k_Target,struct_name);

% Prints .concentration file!
kDiffusion = 1e-6;
print_Concentration_Info(Nx,Nx,Concentration,kDiffusion,struct_name);


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
% FUNCTION: prints MASS points to a file called struct_name.mass
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
function print_Lagrangian_Mass_Pts(xLag,kMass,Mass,struct_name,Nbefore,Nspider)

    %LOOP OVER LAG PTS FOR MASS PTS IN LAGRANGIAN INDEXING

    N = ( length(xLag)-Nbefore ) / Nspider;

    mass_fid = fopen([struct_name '.mass'], 'w');

    fprintf(mass_fid, '%d\n', N );

    %Loops over all Lagrangian Pts.
    %for s = Nbefore+1:length(xLag)
    for s = 1:N
        fprintf(mass_fid, '%d %1.16e %1.16e\n', 1+Nbefore+(s-1)*Nspider, kMass,Mass);
    end

    fclose(mass_fid); 
   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints SPRING points to a file called 'struct_name'.spring
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Springs(xLag,yLag,k_Spring,ds_Rest,struct_name,Nbefore,Nspider)

    N = ( length(xLag)-Nbefore ) / Nspider; % # of total spiders
    Nsprings = N*(Nspider-1);
    
    
    spring_fid = fopen([struct_name '.spring'], 'w');

    fprintf(spring_fid, '%d\n', Nsprings );    % Print # of springs 

    %SPRINGS BETWEEN VERTICES
    for j = 1:N
        for s = 1:(Nspider-1)    
            s1 = 1+Nbefore+(j-1)*Nspider;
            s2 = 2+Nbefore+(j-1)*Nspider;
            fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s1, s2, k_Spring, ds_Rest);  
        end
    end
    fclose(spring_fid);     
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints TARGET points to a file called 'struct_name'.vertex
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
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints BEAM (Torsional Spring) points to a file called 'struct_name'.beam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Beams(xLag,yLag,k_Beam,C,struct_name,Nbefore,Nspider)

    % k_Beam: beam stiffness
    % C: beam curvature
    
    N = ( length(xLag)-Nbefore ) / Nspider % # of total spiders
    Nbeams = N*(Nspider-2)
    
    beam_fid = fopen([struct_name '.beam'], 'w');

    fprintf(beam_fid, '%d\n', Nbeams );

    %spring_force = kappa_spring*ds/(ds^2);

    %BEAMS BETWEEN VERTICES
    for j = 1:N                 %loops over all spiders
        for s = 1:(Nspider-2)    
                s1 = 1+Nbefore+(j-1)*Nspider;
                s2 = 2+Nbefore+(j-1)*Nspider;
                s3 = 3+Nbefore+(j-1)*Nspider;
                %fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s1, s2, k_Spring, ds_Rest);  
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',s1, s2, s3, k_Beam, 0 );  
        end
    end

    fclose(beam_fid); 
   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints CONCENTRATION INFO to file called
%           'struct_name'.concentration
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function print_Concentration_Info(Nx,Ny,C,kDiffusion,struct_name)

    con_fid = fopen([struct_name '.concentration'], 'w');

    fprintf(con_fid, '%d\n', kDiffusion );

    for i=1:Ny
        for j=1:Nx
            fprintf(con_fid, '%1.16e ', C(i,j) );
        end
        fprintf(con_fid,'\n');
    end    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes "curvature" of starting configuration
% 
% NOTE: not curvature in the traditional geometric sense, in the 'discrete'
% sense through cross product.
%
% NOTE: assumes a CLOSED structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C = compute_Curvatures(xLag,yLag)

N = length(xLag);
C = zeros( N );

%Note: needs to be done same order as you print .beam file!
for i=1:N
   
   % Pts Xp -> Xq -> Xr (same as beam force calc.)
   
   if ( (i > 1) && (i < N) )
   
        Xp = xLag(i-1); Xq = xLag(i); Xr = xLag(i+1);
        Yp = yLag(i-1); Yq = yLag(i); Yr = yLag(i+1);
   
   elseif (i==1)
       
        Xp = xLag(N); Xq = xLag(i); Xr = xLag(i+1);
        Yp = yLag(N); Yq = yLag(i); Yr = yLag(i+1);
       
   elseif (i==N)
       
        Xp = xLag(N-1); Xq = xLag(N); Xr = xLag(1);
        Yp = yLag(N-1); Yq = yLag(N); Yr = yLag(1);
       
   end
       
   C(i) = (Xr-Xq)*(Yq-Yp) - (Yr-Yq)*(Xq-Xp); %Cross product btwn vectors
      
end
 





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives initial concentration gradient inside channel
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C = give_Me_Initial_Concentration(Lx,Nx,dx,buffx,buffy)

%WHERE OUTER TUBE LIES
%xMin = 0.15; xMax = 0.45;
%yMin = 0.85; yMax = 1.15;

xMin = buffx; xMax = Lx-buffx;
yMin = buffy; yMax = Lx-buffy;

% xMin = 0; xMax = Lx;
% yMin = 0; yMax = Lx;

xMid = (xMin+xMax)/2;
yMid = (yMin+yMax)/2;

xDiff = (xMax-xMin)/2;
yDiff = (yMax-yMin)/2;

x = 0:dx:Lx;
y = 0:dx:Lx;
inds = give_Me_Indices_To_Apply_Force(x,y,xMin,xMax,yMin,yMax);

%inds(1:length(inds(:,1))/2,:)
%inds(length(inds(:,1))/2+1:end,:)

C = ones(Nx,Nx);
for i=1:length( inds(:,1) )
    xInd = inds(i,1);
    yInd = inds(i,2);
    xPt = x(xInd);
    yPt = y(yInd);
    %C(xInd,yInd ) = (-0.5/yDiff^2)*( (yPt-yMid) - yDiff )*( (yPt-yMid) + yDiff ) +  (-0.5/xDiff^2)*( (xPt-xMid) - xDiff )*( (xPt-xMid) + xDiff ); %1.0;
    %C(yInd,xInd ) = (-1.0/xDiff^2)*( (xPt-xMid) - xDiff )*( (xPt-xMid) + xDiff ); %1.0;
    
    C(yInd,xInd) = 0.5*(tanh(80.0*(yPt-0.5*Lx-0.001*cos(2*pi*xPt)))+1);
    
%     if yPt >= Lx/2+1*dx
%         C(yInd,xInd) = 1;
%     elseif yPt <= Lx/2-1*dx;
%         C(yInd,xInd) = 0;
%     else
%         C(yInd,xInd) = 0.5;
%     end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes indices for placing initial concentration 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function inds = give_Me_Indices_To_Apply_Force(x,y,xMin,xMax,yMin,yMax)

j=1; noMinYet = 1;
while noMinYet
    
    if ( x(j) >= xMin )
        iX_min = j;
        noMinYet = 0;
    end
    j=j+1;
end

j=length(x); noMaxYet = 1;
while noMaxYet
    
    if ( x(j) <= xMax )
        iX_max = j;
        noMaxYet = 0;
    end
    j=j-1;
end

j=1; noMinYet = 1;
while noMinYet
    
    if ( y(j) >= yMin )
        iY_min = j;
        noMinYet = 0;
    end
    j=j+1;
end

j=length(y); noMaxYet = 1;
while noMaxYet
    
    if ( y(j) <= yMax )
        iY_max = j;
        noMaxYet = 0;
    end
    j=j-1;
end

iX_Vec = iX_min:1:iX_max;
iY_Vec = iY_min:1:iY_max;

n = 1;
for i=1:length(iX_Vec)
    for j=1:length(iY_Vec)
        inds(n,1) = iX_Vec(i);
        inds(n,2) = iY_Vec(j);
        n = n+1; 
    end
end





    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates the Lagrangian structure geometry
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [xLag,yLag,xS_all,yD_all,C,Nspider] = give_Me_Immsersed_Boundary_Geometry_and_Concentration(ds,Nx,Lx,dx)
 
% ds: Lagrangian pt. spacing
% Nx: Eulerian grid resolution

Buffx = 0.3*Lx;
Buffy = 0.05*Lx;

fprintf('\nThis is the buffer on x and y, respectively: %d and %d\n',Buffx,Buffy);

ySide = Buffy:ds:Lx-Buffy;
xLeft = Buffx*ones(1,length(ySide));
xRight= (Lx-Buffx)*ones(1,length(ySide));

xBottom = Buffx+ds:ds:Lx-Buffx-ds;
yBottom = Buffy*ones(1,length(xBottom));
yTop = (Lx-Buffy)*ones(1,length(xBottom));

xLag = [xLeft xBottom xRight xBottom];
yLag = [ySide yBottom ySide yTop];

plot(xLeft,ySide,'b*'); hold on;
plot(xBottom,yBottom,'g*'); hold on;
plot(xRight,ySide,'r*'); hold on;
axis([0 Lx 0 Lx]); 
 
%xSpider = [0.25 0.3 0.35 0.4 0.45 0.48 0.55 0.6 0.65 0.7 0.75];
% x-spider coords
xS = 0.325:0.05:0.675;
%yS = 0.5*ones(1,length(xSpider));

% Original Drag Line
yD = 0.5:ds:0.625;

% Length in each spider
Nspider = length(yD);

yD_all = [];
for i=1:length(xS)
   yD_all = [yD_all yD]; 
end

xS_all = [];
for j=1:length(xS)
    x = xS(j);
    for i=1:Nspider
        xS_all = [xS_all x];
    end
end

plot(xS_all,yD_all,'k*'); hold on;
%plot(xSpider,ySpider,'k*'); hold on;

C = give_Me_Initial_Concentration(Lx,Nx,dx,Buffx,Buffy);





