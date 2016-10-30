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
% If you would like us %to add a specific muscle model, please let Nick (nick.battista@unc.edu) know.
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates the BABY_SPIDER-EXAMPLE geometry and prints associated input files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Baby_Spider()

%
% Grid Parameters (MAKE SURE MATCHES IN input2d !!!)
%
Nx =  256;        % # of Eulerian Grid Pts. in x-Direction (MUST BE EVEN!!!)
Ny =  256;        % # of Eulerian Grid Pts. in y-Direction (MUST BE EVEN!!!)
Lx = 3.0;        % Length of Eulerian Grid in x-Direction
Ly = 3.0;        % Length of Eulerian Grid in y-Direction
dx = Lx/Nx;      % Grid Resolution


% Immersed Structure Geometric / Dynamic Parameters %
Ls = 1.0;        % Length of baby-spider web
ds = Lx/(2*Nx);  % Lagrangian spacing, ds
struct_name = 'spider'; % Name for .vertex, .spring, etc files.


% Call function to construct geometry (spider)
[xLagSpider,yLagSpider,Ninfo] = give_Me_Immsersed_Boundary_Geometry(ds,Ny,Ly,Ls);

% Call function to construct geometry (outer square)
[xLagSq,yLagSq,Concentration] = give_Me_Immsersed_Boundary_Geometry_and_Concentration(ds,Nx,Lx,dx);


% Plot Geometry to test
plot(xLagSpider,yLagSpider,'r-'); hold on;
plot(xLagSpider,yLagSpider,'*'); hold on;
%plot(xLagSq,yLagSq,'k'); hold on;
xlabel('x'); ylabel('y');
axis([0 Lx 0 Ly]);


% Print INFO to screen
fprintf('\n\n                   INFO TO CHECK INPUT FILES!\n\n');
%fprintf('# of Lag. Pts in WEB: %d\n',Ninfo(1));
%fprintf('# of Lag. Pts in FLOOR: %d\n',Ninfo(2));
%fprintf('Total # of Lag. Pts: %d\n',Ninfo(1)+Ninfo(2));
fprintf('Index of MASSIVE Pt (last pt. of WEB): %d\n\n\n',Ninfo(1));


% Combine Geometries
%offset = length(xLagSpider);
%xLag = [xLagSpider xLagSq];
%yLag = [yLagSpider yLagSq];
xLag = xLagSpider;
yLag = yLagSpider;

% Prints .vertex file!
print_Lagrangian_Vertices(xLag,yLag,struct_name);


% Prints .mass file!
k_Mass = 6.82e4;        % 'spring' stiffness parameter for tethering
Mass =   2.0e-3;        % "MASS" value for 'MASSIVE' nodal movement
print_Lagrangian_Mass_Pts(xLagSpider,k_Mass,Mass,struct_name,Ninfo);


% Prints .spring file!
k_Spring = 6.82e4;
print_Lagrangian_Springs(xLagSpider,yLag,k_Spring,struct_name,Ninfo,ds);


% Prints .beam file!
k_Beam = 1.98259e3; C = 0.0;
print_Lagrangian_Beams(xLagSpider,yLagSpider,k_Beam,C,struct_name,Ninfo);


% Prints .concentration file!
kDiffusion = 1e-6;
print_Concentration_Info(Nx,Nx,Concentration,kDiffusion,struct_name);


% Prints .target file!
%k_Target = 1e7;
%print_Lagrangian_Target_Pts(xLagSq,k_Target,struct_name,offset);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints VERTEX points to a file called spider.vertex
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
% FUNCTION: prints MASS points to a file called spider.mass
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
function print_Lagrangian_Mass_Pts(xLag,kMass,Mass,struct_name,Ninfo)

    %LOOP OVER LAG PTS FOR MASS PTS IN LAGRANGIAN INDEXING

    N = 1;  % ONLY 1 MASS PT.!

    mass_fid = fopen([struct_name '.mass'], 'w');

    fprintf(mass_fid, '%d\n', N );

    % Single Mass Pt. (bottom point of web!)
    fprintf(mass_fid, '%d %1.16e %1.16e\n', Ninfo(1), kMass,Mass);
    
    fclose(mass_fid); 
       
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints Vertex points to a file called spider.vertex
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Target_Pts(xLag,k_Target,struct_name,Ninfo,offset)

    N = length(xLag);
    
    target_fid = fopen([struct_name '.target'], 'w');

    fprintf(target_fid, '%d\n', N );

    %Loops over Lag. Pts on FLOOR ONLY
    for s = 1:N
        fprintf(target_fid, '%d %1.16e\n', s+offset, k_Target);
    end

    fclose(target_fid); 
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints BEAM (Torsional Spring) points to a file called spider.beam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Beams(xLag,yLag,k_Beam,C,struct_name,Ninfo)

    % k_Beam: beam stiffness
    % C: beam curvature
    
    N = Ninfo(1)-2; % NOTE: Total # of beams = % of Lag Pts. on WEB - 2

    beam_fid = fopen([struct_name '.beam'], 'w');

    fprintf(beam_fid, '%d\n', N );

    %BEAMS BETWEEN VERTICES
    for s = 2:N-1
        fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',s-1, s, s+1, k_Beam, C);  
    end
    fclose(beam_fid); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints SPRING points to a file called spider.spring
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Springs(xLag,yLag,k_Spring,struct_name,Ninfo,ds)

    N = Ninfo(1)-1;  % # of springs (just on 'web' of geometry)

    spring_fid = fopen([struct_name '.spring'], 'w');

    fprintf(spring_fid, '%d\n', N );

    %spring_force = kappa_spring*ds/(ds^2);

    %SPRINGS BETWEEN VERTICES ON WEBBING ONLY
    for s = 1:N
        fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, s+1, k_Spring, ds);  
    end
    fclose(spring_fid); 
    
    
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
% FUNCTION: creates the Lagrangian structure geometry
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag,Ninfo] = give_Me_Immsersed_Boundary_Geometry(ds,Ny,Ly,Ls)

% ASSUMES Nx=Ny %

%
% Create Spider Web Geometry (springs, beams, and 1-mass pt)
%
yS =   2.0;                       % Highest Point of Web
yWeb = yS:-ds:yS-Ls;                  % yPts of web
Nweb = length(yWeb);                  % # of Lag. Pts in Web
xWeb = (1.25)*ones(1,Nweb);           % xPts of web

xLag = xWeb; yLag = yWeb;
Ninfo(1) = Nweb;

%
% Create Floor (for target pts)
% 
% frac = 0.2;                                       % Fraction of bottom of domain, 'not floored'
% frac2 = 0.05;                                     % Fraction of height of domain for floor to be placed
% xFloor = (frac/2)*Ly:ds:Ly*( 1.0 - (frac/2) );    % xPts of floor
% Nfloor = length(xFloor);                          % # of Lag. Pts in Floor
% yFloor = frac2*Ly*ones(1,Nfloor);                 % yPts of floor
% 
% 
% xLag = [xWeb xFloor];             %Vector of all x-Lag. Pts.
% yLag = [yWeb yFloor];             %Vector of all y-Lag. Pts.
% Ninfo= [Nweb Nfloor Nweb+Nfloor]; %Vector of # of Lag. Pts for each part of geometry



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives initial concentration gradient inside channel
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C = give_Me_Initial_Concentration(Lx,Nx,dx,buffx,buffy)

%WHERE OUTER TUBE LIES

% xMin = buffx; xMax = Lx-buffx;
% yMin = buffy; yMax = Lx-buffy;

xMin = 0; xMax = Lx;
yMin = 0; yMax = Lx;

xMid = (xMin+xMax)/2;
yMid = (yMin+yMax)/2;

xDiff = (xMax-xMin)/2;
yDiff = (yMax-yMin)/2;

x = 0:dx:Lx;
y = 0:dx:Lx;
inds = give_Me_Indices_To_Apply_Force(x,y,xMin,xMax,yMin,yMax);

C = ones(Nx,Nx);
for i=1:length( inds(:,1) )
    xInd = inds(i,1);
    yInd = inds(i,2);
    xPt = x(xInd);
    yPt = y(yInd);
    C(yInd,xInd) = 0.02*yPt;
    %C(xInd,yInd ) = (-0.5/yDiff^2)*( (yPt-yMid) - yDiff )*( (yPt-yMid) + yDiff ) +  (-0.5/xDiff^2)*( (xPt-xMid) - xDiff )*( (xPt-xMid) + xDiff ); %1.0;
    %C(yInd,xInd ) = (-1.0/xDiff^2)*( (xPt-xMid) - xDiff )*( (xPt-xMid) + xDiff ); %1.0;
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
 
function [xLag,yLag,C] = give_Me_Immsersed_Boundary_Geometry_and_Concentration(ds,Nx,Lx,dx)
 
% ds: Lagrangian pt. spacing
% Nx: Eulerian grid resolution

Buffx = 0.025*Lx;
Buffy = 0.025*Lx;

%fprintf('\nThis is the buffer on x and y, respectively: %d and %d\n',Buffx,Buffy);

ySide = Buffy:ds:Lx-Buffy;
xLeft = Buffx*ones(1,length(ySide));
xRight= (Lx-Buffx)*ones(1,length(ySide));

xBottom = Buffx+ds:ds:Lx-Buffx-ds;
yBottom = Buffy*ones(1,length(xBottom));
yTop = (Lx-Buffy)*ones(1,length(xBottom));

xLag = [xLeft xBottom xRight xBottom];
yLag = [ySide yBottom ySide yTop];

%plot(xLeft,ySide,'b*'); hold on;
%plot(xBottom,yBottom,'g*'); hold on;
%plot(xRight,ySide,'r*'); hold on;
%axis([0 Lx 0 Lx]); 
 
C = give_Me_Initial_Concentration(Lx,Nx,dx,Buffx,Buffy);

