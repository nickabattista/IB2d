%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  battistn@tcnj.edu
% Date Created: September 9th, 2016
% Institution Created: UNC-CH
% Date Modified: June 26, 2021
% Institution: TCNJ
%
% This code is capable of creating Lagrangian Structures using:
% 	1. Springs
% 	2. Beams (*torsional springs*)
% 	3. Target Points
%	4. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-Tension)")
%
% One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting lengths, etc
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
% 
% If you would like us to add a specific muscle model, please let Nick (battistn@tcnj.edu) know.
%
% If you use this code for the purposes of teaching, research, or recreation please let Nick know as well :)
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
Nx = 1024;       % # of Eulerian Grid Pts. in x-Direction (MUST BE EVEN!!!)
Ny = 512;        % # of Eulerian Grid Pts. in y-Direction (MUST BE EVEN!!!)
Lx = 1.0;        % Length of Eulerian Grid in x-Direction
dx = Lx/Nx;      % Grid spatial resolution

%
% Immersed Structure Geometric / Dynamic Parameters %
%
ds = 0.25*dx;              % Lagrangian Pt. Spacing (2x resolution of Eulerian grid)
struct_name = 'boussinesq'; % Name for .vertex, .spring, etc files. (must match what's in 'input2d')


% Call function to construct geometry
[xLag,yLag,Concentration,inds,OffY] = give_Me_Immsersed_Boundary_Geometry_and_Concentration(ds,Nx,Ny,Lx,dx);

size(Concentration)


xLag = dx/2:dx/2:Lx-dx/2;
yLag = 0.01*ones(1,length(xLag));
yLag2 = 0.49*ones(1,length(xLag));

xLag = [xLag xLag];
yLag = [yLag yLag2];

%--------------------------------------
% Lines of Tracer Particles
%--------------------------------------
OffY=0;
xLagT = 0:2*ds:Lx;
yLag1 = (0.15+OffY)*ones( size( xLagT ) );
yLag2 = (0.25+OffY)*ones( size( xLagT ) );
yLag3 = (0.35+OffY)*ones( size( xLagT ) );
%

%--------------------------------------
% Line at Bottom of Domain
%--------------------------------------
xLagB = 0:0.75*ds:Lx;
yLagB1 = 1/3*ds*ones( size( xLagB ) );
yLagB2 = 2/3*ds*ones( size( xLagB ) );
yLagB3 = ds*ones( size( xLagB ) );

%--------------------------------------
% Combine Geometries
%--------------------------------------
%xLag = [xLagB xLagB xLagB];
%yLag = [yLagB1 yLagB2 yLagB3];
%
xLag = [xLagT xLagT xLagT];
yLag = [yLag1 yLag2 yLag3];

% Plot Geometry to test
plot(xLag,yLag,'r-'); hold on;
plot(xLag,yLag,'*'); hold on;
xlabel('x'); ylabel('y');
axis([0 Lx 0 Lx]); 



% Prints .vertex file!
print_Lagrangian_Vertices(xLag,yLag,struct_name);


% Prints .spring file!
%k_Spring = 2.5e4;                    % Spring stiffness (does not need to be equal for all springs)
%ds_Rest = 0.0;                       % Spring resting length (does not need to be equal for all springs)
%print_Lagrangian_Springs(xLag,yLag,k_Spring,ds_Rest,struct_name);


% Prints .beam file!
% k_Beam = 0.5;                      % Beam Stiffness (does not need to be equal for all beams)
% C = compute_Curvatures(xLag,yLag)  % Computes curvature of initial configuration
%print_Lagrangian_Beams(xLag,yLag,k_Beam,C,struct_name);


% Prints .target file!
k_Target = 1e6; %1e7
print_Lagrangian_Target_Pts(xLag,k_Target,struct_name);

% Prints .concentration file!
print_Concentration_Info(Nx,Ny,Concentration,struct_name);

% Print .mesh file where Boussinesq is on Eulerian Mesh
print_Boussinesq_Mesh(inds)

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
% FUNCTION: prints BOUSSINESQ affected EULERIAN point indices
%           boussinesq.mesh
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Boussinesq_Mesh(inds)

    N = length(inds(:,1)); % Total # of Affected Eulerian Pts

    bouss_fid = fopen('boussinesq.mesh', 'w');

    %Loops over all Lagrangian Pts.
    for s = 1:N
        xInd = inds(s,1);
        yInd = inds(s,2);
        
        if ( ( yInd < 0.425 ) || ( yInd > 0.025 ) )
            fprintf(bouss_fid, '%d %d\n', xInd,yInd);
        end
    end

    fclose(bouss_fid);
   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints SPRING points to a file called 'struct_name'.spring
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Springs(xLag,yLag,k_Spring,ds_Rest,struct_name)

    N = length(xLag);

    spring_fid = fopen([struct_name '.spring'], 'w');

    fprintf(spring_fid, '%d\n', N );    % Print # of springs 

    %SPRINGS BETWEEN VERTICES
    for s = 1:N
            if s < N         
                fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, s+1, k_Spring, ds_Rest);  
            else
                %Case s=N
                fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, 1,   k_Spring, ds_Rest);  
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

function print_Lagrangian_Beams(xLag,yLag,k_Beam,C,struct_name)

    % k_Beam: beam stiffness
    % C: beam curvature
    
    N = length(xLag); % NOTE: Total number of beams = Number of Total Lag Pts. - 2

    beam_fid = fopen([struct_name '.beam'], 'w');

    fprintf(beam_fid, '%d\n', N );

    %spring_force = kappa_spring*ds/(ds^2);

    %BEAMS BETWEEN VERTICES
    for s = 2:N-1
            if  s <= N-1         
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',s-1, s, s+1, k_Beam, C(s) );  
            else
                %Case s=N
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',s-1, s, 1,   k_Beam, C(s) );  
            end
    end
    fclose(beam_fid); 
   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints CONCENTRATION INFO to file called
%           'struct_name'.concentration
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function print_Concentration_Info(Nx,Ny,C,struct_name)

    con_fid = fopen([struct_name '.concentration'], 'w');


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

function [C,inds,OffY] = give_Me_Initial_Concentration(Lx,Nx,Ny,dx,buffx,buffy)


%xMin = buffx; xMax = Lx-buffx;
xMin = 0; xMax = Lx;

yMin = 0; yMax = Lx/2;



x=linspace(0,Lx,Nx);
y=linspace(0,Lx/2,Ny);
inds = give_Me_Indices_To_Apply_Force(x,y,xMin,xMax,yMin,yMax);

%inds(1:length(inds(:,1))/2,:)
%inds(length(inds(:,1))/2+1:end,:)

C = ones(Ny,Nx);
for i=1:length( inds(:,1) )
    xInd = inds(i,1);
    yInd = inds(i,2);
    xPt = x(xInd);
    yPt = y(yInd);
    %C(xInd,yInd ) = (-0.5/yDiff^2)*( (yPt-yMid) - yDiff )*( (yPt-yMid) + yDiff ) +  (-0.5/xDiff^2)*( (xPt-xMid) - xDiff )*( (xPt-xMid) + xDiff ); %1.0;
    %C(yInd,xInd ) = (-1.0/xDiff^2)*( (xPt-xMid) - xDiff )*( (xPt-xMid) + xDiff ); %1.0;
    
    %C(yInd,xInd) = 0.5*(tanh(80.0*(yPt-0.25*Lx-0.001*cos(pi*xPt)))+1);
%     if yPt > (0.25+0.025*sin(20*pi*xPt))
%         C(yInd,xInd) = 0;
%     else
%         C(yInd,xInd) = 1;
%     end

    %-----------------------------------------
    % Makes 4 equal-width horizontal strips
    %-----------------------------------------
    Scale = 1;
    Amp = 0.025;
    OffY = Amp;
    lambda = 0.2;
    %
    if yPt < ( 0.125+Amp*sin(2*pi*xPt/lambda) + OffY ) 
        C(yInd,xInd) = 0;
    elseif yPt < ( 0.25+Amp*sin(2*pi*xPt/lambda) + 0 )
        C(yInd,xInd) = 1*Scale;
    elseif yPt < ( 0.375+Amp*sin(2*pi*xPt/lambda) - OffY )
        C(yInd,xInd) = 2.0*Scale;
    else
        C(yInd,xInd) = 3.0*Scale;
    end
    

    %-----------------------------------------
    % Makes 4 Diagonally-aligned strips
    %-----------------------------------------
%     % Define Max Vertical Values on interval
%     y0 = 0;
%     y1 = 0.125;
%     y2 = 0.25;
%     y3 = 0.375;
%     y4 = 0.5;
%     
%     % Define Slopes
%     m4 = (y3-y4);
%     m3 = (y2-y3);
%     m2 = (y1-y2);
%     m1 = (y0-y1);
%     
%     % Define y value along lines at specific xPt
%     yL_1 = m1*( xPt ) + y1;
%     yL_2 = m2*( xPt ) + y2;
%     yL_3 = m3*( xPt ) + y3;
%     yL_4 = m4*( xPt ) + y4;
%     
%     if yPt < yL_1
%         C(yInd,xInd) = 0;
%     elseif yPt < yL_2
%         C(yInd,xInd) = 1;
%     elseif yPt < yL_3
%         C(yInd,xInd) = 2;
%     elseif yPt < yL_4
%         C(yInd,xInd) = 3;
%     else
%         C(yInd,xInd) = 4;
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
 
function [xLag,yLag,C,inds,OffY] = give_Me_Immsersed_Boundary_Geometry_and_Concentration(ds,Nx,Ny,Lx,dx)
 
% ds: Lagrangian pt. spacing
% Nx: Eulerian grid resolution

Buffx = 0;%0.05*Lx;
Buffy = 0;%0.05*Lx;

fprintf('\nThis is the buffer on x and y, respectively: %d and %d\n',Buffx,Buffy);

ySide = Buffy:ds:Lx/2-Buffy;
xLeft = Buffx*ones(1,length(ySide));
xRight= (Lx-Buffx)*ones(1,length(ySide));

xBottom = Buffx+ds:ds:Lx-Buffx-ds;
yBottom = Buffy*ones(1,length(xBottom));
yTop = (Lx/2-Buffy)*ones(1,length(xBottom));

xLag = [xLeft xBottom xRight xBottom];
yLag = [ySide yBottom ySide yTop];

% plot(xLeft,ySide,'b*'); hold on;
% plot(xBottom,yBottom,'g*'); hold on;
% plot(xRight,ySide,'r*'); hold on;
% axis([0 Lx 0 Lx]); 
 
[C,inds,OffY] = give_Me_Initial_Concentration(Lx,Nx,Ny,dx,Buffx,Buffy);





