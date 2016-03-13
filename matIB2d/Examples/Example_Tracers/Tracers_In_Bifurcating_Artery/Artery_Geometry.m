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
% FUNCTION: creates the BIFURCATING_ARTERY-EXAMPLE geometry and prints associated input files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Artery_Geometry()

%
% Grid Parameters (MAKE SURE MATCHES IN input2d !!!)
%
Nx =  128;        % # of Eulerian Grid Pts. in x-Direction (MUST BE EVEN!!!)
Ny =  128;        % # of Eulerian Grid Pts. in y-Direction (MUST BE EVEN!!!)
Lx = 2.0;        % Length of Eulerian Grid in x-Direction
Ly = 2.0;        % Length of Eulerian Grid in y-Direction


% Immersed Structure Geometric / Dynamic Parameters %
ds= min(Lx/(2*Nx),Ly/(2*Ny));  % Lagrangian spacing
L = 0.25*Lx;                    % Length of Channel
w = 0.15*Ly;                    % Width of Channel
x0 = 0.4;                      % x-Center for Cylinder
y0 = 1.0;                      % y-Center for Cylinder
r = w/6;                       % Radii of Cylinder
struct_name = 'artery'; % Name for .vertex, .spring, etc files.


% Call function to construct BASE CHANNEL geometry
[xLag,yLag] = give_Me_Channel_Immsersed_Boundary_Geometry(ds,L,w,Lx,Ly);
[xLag_C,yLag_C] = give_Me_Cylinder_Immsersed_Boundary_Geometry(ds,r,x0,y0);

% Call functions to construct FIRST bifurcation
xlast = xLag(end);  ylast_top = yLag(end); ylast_bot = yLag(end/2);
h = Ly/10;
w1 = 0.425*w;
w2 = 0.425*w;
[xBif1,yBif1] = compute_Bifurcation_Tanh(ds,h,L,xlast,ylast_top,ylast_bot,w1,w2);
dia = yBif1(3*end/4+1) - yBif1(end/4+1); 
yMid = ( yBif1(3*end/4+1) + yBif1(end/4+1) ) / 2;
xMid = xBif1(3*end/4+1);
[xR1,yR1] = compute_Bifurcation_Cap(3*ds/4,dia,xMid,yMid);


% Call functions to construct TOP-RIGHT bifurcation
xlast = xBif1(end);  ylast_top = yBif1(3*end/4); ylast_bot = yBif1(end);
h = Ly/15;
w1 = 0.425*(0.425*w);
w2 = 0.425*(0.425*w);
[xBif2,yBif2] = compute_Bifurcation_Tanh(ds,h,L,xlast,ylast_top,ylast_bot,w1,w2);
dia = yBif2(3*end/4+1) - yBif2(end/4+1); 
yMid = ( yBif2(3*end/4+1) + yBif2(end/4+1) ) / 2;
xMid = xBif2(3*end/4+1);
[xR2,yR2] = compute_Bifurcation_Cap(ds/2,dia,xMid,yMid);

% Call functions to construct BOTTOM-RIGHT bifurcation
xlast = xBif1(end/2);  ylast_top = yBif1(end/2); ylast_bot = yBif1(end/4);
h = Ly/15;
w1 = 0.425*(0.425*w);
w2 = 0.425*(0.425*w);
[xBif3,yBif3] = compute_Bifurcation_Tanh(ds,h,L,xlast,ylast_top,ylast_bot,w1,w2);
dia = yBif3(3*end/4+1) - yBif3(end/4+1); 
yMid = ( yBif3(3*end/4+1) + yBif3(end/4+1) ) / 2;
xMid = xBif3(3*end/4+1);
[xR3,yR3] = compute_Bifurcation_Cap(ds/2,dia,xMid,yMid);


% Call function for placing tracer points
[xT,yT] = give_Me_Tracer_Points(ds,r,x0,y0,w,L);


% Plot Geometry to test
plot(xLag(1:end/2),yLag(1:end/2),'r-'); hold on;
plot(xLag(end/2+1:end),yLag(end/2+1:end),'r-'); hold on;
%
plot(xBif1(1:end/4),yBif1(1:end/4),'k-'); hold on;
plot(xBif1(end/4+1:end/2),yBif1(end/4+1:end/2),'k-'); hold on;
plot(xBif1(end/2+1:3*end/4),yBif1(end/2+1:3*end/4),'k-'); hold on;
plot(xBif1(3*end/4+1:end),yBif1(3*end/4+1:end),'k-'); hold on;
%
plot(xR1,yR1,'m*'); hold on;
plot(xR2,yR2,'m*'); hold on;
plot(xR3,yR3,'m*'); hold on;
%
plot(xBif2(1:end/4),yBif2(1:end/4),'k-'); hold on;
plot(xBif2(end/4+1:end/2),yBif2(end/4+1:end/2),'k-'); hold on;
plot(xBif2(end/2+1:3*end/4),yBif2(end/2+1:3*end/4),'k-'); hold on;
plot(xBif2(3*end/4+1:end),yBif2(3*end/4+1:end),'k-'); hold on;
%
plot(xBif3(1:end/4),yBif3(1:end/4),'k-'); hold on;
plot(xBif3(end/4+1:end/2),yBif3(end/4+1:end/2),'k-'); hold on;
plot(xBif3(end/2+1:3*end/4),yBif3(end/2+1:3*end/4),'k-'); hold on;
plot(xBif3(3*end/4+1:end),yBif3(3*end/4+1:end),'k-'); hold on;
%
%plot(xLag_C,yLag_C,'r-'); hold on;

plot(xLag,yLag,'*'); hold on;
%plot(xLag_C,yLag_C,'g*'); hold on;
plot(xBif1,yBif1,'r*'); hold on;
plot(xT,yT,'m*'); hold on;
xlabel('x'); ylabel('y');
axis([0 Lx 0 Ly]);




%Bottom Branches
XBot = [xLag(1:end/2) xBif1(1:end/4) xBif3(1:end/4) xBif3(end/2:-1:end/4+1) xR3(end:-1:1) xBif3(3*end/4+1:end) xBif3(3*end/4:-1:end/2+1) xBif1(end/2:-1:end/4+1) xR1(end:-1:1) ];
YBot = [yLag(1:end/2) yBif1(1:end/4) yBif3(1:end/4) yBif3(end/2:-1:end/4+1) yR3(end:-1:1) yBif3(3*end/4+1:end) yBif3(3*end/4:-1:end/2+1) yBif1(end/2:-1:end/4+1) yR1(end:-1:1) ];

XTop = [xBif1(3*end/4+1:end) xBif2(1:end/4) xBif2(end/2:-1:end/4+1) xR2(end:-1:1) xBif2(3*end/4+1:end) xBif2(3*end/4:-1:end/2+1) xBif1(3*end/4:-1:end/2+1) xLag(end:-1:end/2+1) ];
YTop = [yBif1(3*end/4+1:end) yBif2(1:end/4) yBif2(end/2:-1:end/4+1) yR2(end:-1:1) yBif2(3*end/4+1:end) yBif2(3*end/4:-1:end/2+1) yBif1(3*end/4:-1:end/2+1) yLag(end:-1:end/2+1) ];

xLag = [XBot XTop];
yLag = [YBot YTop];

% Test Ordering!
figure(2)
for i=1:length(xLag)
   plot(xLag(i),yLag(i),'*'); hold on;
   axis([0 Lx 0 Ly]);
   pause(0.0001)
end


% Prints .vertex file!
print_Lagrangian_Vertices(xLag,yLag,struct_name);

% Prints .tracer file!
print_Lagrangian_Tracers(xT,yT,struct_name)

% Prints .spring file!
%k_Spring = 1e7;
%print_Lagrangian_Springs(xLag,yLag,k_Spring,ds_Rest,struct_name);


% Prints .beam file!
%k_Beam = 0.5; C = 0.0;
%print_Lagrangian_Beams(xLag,yLag,k_Beam,C,struct_name);


% Prints .target file!
k_Target = 1.25e7;
print_Lagrangian_Target_Pts(xLag,k_Target,struct_name);

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
% FUNCTION: prints TRACER points to a file called rubberband.vertex
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Tracers(xLag,yLag,struct_name)

    N = length(xLag);

    vertex_fid = fopen([struct_name '.tracer'], 'w');

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
% FUNCTION: prints TARGET points to a file called rubberband.vertex
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
% FUNCTION: prints BEAM (Torsional Spring) points to a file called rubberband.beam
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
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',s-1, s, s+1, k_Beam, C);  
            else
                %Case s=N
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e\n',s-1, s, 1,   k_Beam, C);  
            end
    end
    fclose(beam_fid); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints SPRING points to a file called rubberband.spring
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Springs(xLag,yLag,k_Spring,ds_Rest,struct_name)

    N = length(xLag);

    spring_fid = fopen([struct_name '.spring'], 'w');

    fprintf(spring_fid, '%d\n', N );

    %spring_force = kappa_spring*ds/(ds^2);

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
% FUNCTION: creates the Lagrangian structure geometry for channel
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag] = give_Me_Channel_Immsersed_Boundary_Geometry(ds,L,w,Lx,Ly)

% The immsersed structure is a channel %
x = (Lx/2-1.7*L)/2:ds:L+(Lx/2-1.7*L)/2;  %xPts
yBot = (Ly-w)/2;               %yVal for bottom of Channel
yTop = Ly - (Ly-w)/2;          %yVal for top of Channel

xLag = [x x];
yLag = [yBot*ones(1,length(x)) yTop*ones(1,length(x))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates the Lagrangian structure geometry for cylinder
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag] = give_Me_Cylinder_Immsersed_Boundary_Geometry(ds,r,x0,y0)

% The immsersed structure is a cylinder %

dtheta = ds/ (2*r);
theta = 0; i=1;
while theta < 2*pi
   xLag(i) = x0 - r*cos(theta);
   yLag(i) = y0 - r*sin(theta);
   theta = theta + dtheta;
   i=i+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives outer bifurcation piece using tanh functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xBif,yBif] = compute_Bifurcation_Tanh(ds,A,L,x0,y0_t,y0_b,w1,w2)

% Compute one leg of branch using tanh function
[x,y] = compute_Tanh_Branch(A,L,ds);

xBif_T1 = L/2 + x0 + ds + x;

yBif_T1 = A + y0_t + y; 
yBif_T2 = A + y0_t + y - w1;

yBif_B1 = y0_b - y - A ;
yBif_B2 = y0_b - y - A + w2;

%plot(xBif_T1,yBif_T1,'g*'); hold on;
%plot(xBif_T1,yBif_B1,'g*'); hold on;
%plot(xBif_T1,yBif_T2,'k*'); hold on;
%plot(xBif_T1,yBif_B2,'k*'); hold on;

xBif = [xBif_T1 xBif_T1 xBif_T1 xBif_T1];
yBif = [yBif_B1 yBif_B2 yBif_T1 yBif_T2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives actually tanh pieces
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,y] = compute_Tanh_Branch(A,L,ds)

%initiate
ct = 1;
xN = -L/2;             x(ct) = xN;
yN = A*tanh(10*xN);    y(ct) = yN;

while xN <= L/2 
    
    %update counter
    ct = ct+1;
    
    xP = x(ct-1);           % x-Prev
    yP = y(ct-1);           % y-Prev
    
    xN = xP;                % left  bound 
    xF = xN+2*ds;           % far guess (right bound)
    xN1 = xN+ds;            % x-guess 
    yN1 = A*tanh(10*xN1);   % y-guess
    err = ( ds - sqrt( (xN1-xN)^2 + (yN1-yN)^2 ) );
    
    while abs(err) > 1e-6
       
       if err > 0
          
          xN = xN1;             % Update LEFT PT. [xN,xN1,xF]
          xN1 = (xN1+xF)/2;     % New Guess
          yN1 = A*tanh(10*xN1); % New Guess
          
          
       elseif err < 0
           
          xF = xN1;             % Update RIGHT PT. [xN,xN1,xF] 
          xN1 = (xN1+xN)/2;     % New Guess
          yN1 = A*tanh(10*xN1); % New Guess
           
       end
        
       %compute error
       err = ( ds - sqrt( (xN1-xP)^2 + (yN1-yP)^2 ) );
        
    end

    %save values
    x(ct) = xN1;
    y(ct) = yN1;
    
    %redefine parameters
    xN = xN1; yN=yN1;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives bifurcation cap, aka half circle.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,y] = compute_Bifurcation_Cap(ds,dia,xMid,yMid)

t = pi/2;
r = dia/2;
dt = ds/r;
n=1;
while t < 3*pi/2
   x(n) = xMid + r*cos(t);
   y(n) = yMid + r*sin(t);
   t = t+dt;
   n = n+1; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives Tracer (x,y) locations (arbitrarily set)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xT,yT] = give_Me_Tracer_Points(ds,r,x0,y0,w,L)

xR = x0 - r -ds - 0.05;
xRM= xR-2*ds;
xM = xR-4*ds;
xLM= xR-6*ds;
xL = xR-8*ds;
xR1= xR+2*ds;
xR2= xR+4*ds;
xR3= xR+6*ds;
xR4= xR+8*ds;
xR5= xR+10*ds;
xR6= xR+12*ds;
xR7= xR+14*ds;
xR8= xR+16*ds;
xR9= xR+18*ds;
xR10= xR+20*ds;
xR11= xR+22*ds;
xR12= xR+24*ds;
xR13= xR+26*ds;
xR14= xR+28*ds;
xR15= xR+30*ds;
xR16= xR+32*ds;

y = (y0-w/2+2*ds:2*ds:y0+w/2-2*ds);

xR = xR*ones(1,length(y));
xRM= xRM*ones(1,length(y));
xM = xM*ones(1,length(y));
xLM= xLM*ones(1,length(y));
xL = xL*ones(1,length(y));
xR1 = xR1*ones(1,length(y));
xR2 = xR2*ones(1,length(y));
xR3 = xR3*ones(1,length(y));
xR4 = xR4*ones(1,length(y));
xR5 = xR5*ones(1,length(y));
xR6 = xR6*ones(1,length(y));
xR7 = xR7*ones(1,length(y));
xR8 = xR8*ones(1,length(y));
xR9 = xR9*ones(1,length(y));
xR10 = xR10*ones(1,length(y));
xR11 = xR11*ones(1,length(y));
xR12 = xR12*ones(1,length(y));
xR13 = xR13*ones(1,length(y));
xR14 = xR14*ones(1,length(y));
xR15 = xR15*ones(1,length(y));
xR16 = xR16*ones(1,length(y));


xT = [xR xRM xM xLM xL xR1 xR2 xR3 xR4 xR5 xR6 xR7 xR8 xR9 xR10 xR11 xR12 xR13 xR14 xR15 xR16];
yT = [y y y y y y y y y y y y y y y y y y y y y];

   
