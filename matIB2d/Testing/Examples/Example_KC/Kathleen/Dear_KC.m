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
% FUNCTION: creates the RUBBERBAND-EXAMPLE geometry and prints associated input files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Dear_KC()

%
% Grid Parameters (MAKE SURE MATCHES IN input2d !!!)
%
Nx =  32;        % # of Eulerian Grid Pts. in x-Direction (MUST BE EVEN!!!)
Ny =  32;        % # of Eulerian Grid Pts. in y-Direction (MUST BE EVEN!!!)
Lx = 1.0;        % Length of Eulerian Grid in x-Direction
Ly = 1.0;        % Length of Eulerian Grid in y-Direction
ds = Lx/(2*Nx);  % Spatial step!


% Immersed Structure Geometric / Dynamic Parameters %
N = 2*Nx;        % Number of Lagrangian Pts. (2x resolution of Eulerian grid)
ds_Rest = 0;     % Resting length of springs
struct_name = 'dear_KC'; % Name for .vertex, .spring, etc files.

% Line on top and bottom of geometry
ds2 = 0.5*ds;
xLine = Lx/20:ds2:0.95*Lx;
yLineB = 0.025*Lx*ones(1,length(xLine));
yLineT= 0.975*Lx*ones(1,length(xLine));

xSq = [xLine xLine yLineB yLineT];
ySq = [yLineB yLineT xLine xLine];

% Call function to construct geometry
[x1,y1] = give_Me_Immsersed_Boundary_Geometry_1(Lx,Nx,ds);
x1 = [x1 xSq];
y1 = [y1 ySq];

% Call function to construct geometry
[x2,y2] = give_Me_Immsersed_Boundary_Geometry_2(Lx,Nx,ds);
x2 = [x2 xSq];
y2 = [y2 ySq];

% Call function to construct geometry
[x3,y3] = give_Me_Immsersed_Boundary_Geometry_3(Lx,Nx,ds);
x3 = [x3 xSq];
y3 = [y3 ySq];

% Plot Geometry to test BEFORE taking out pts.
figure(1)
plot(x1,y1,'r*'); hold on;
axis([0 Lx 0 Ly]);

figure(2)
plot(x2,y2,'b*'); hold on;
axis([0 Lx 0 Ly]);

figure(3)
plot(x3,y3,'k*'); hold on;
axis([0 Lx 0 Ly]);


% Print files to .txt files
please_Print_Vertices_To_File(x1,y1,x2,y2,x3,y3)

% Prints .vertex file!
print_Lagrangian_Vertices(x1,y1,struct_name);


% Prints .spring file!
%k_Spring = 1e7;
%print_Lagrangian_Springs(xLag,yLag,k_Spring,ds_Rest,struct_name);


% Prints .beam file!
%k_Beam = 0.5; C = 0.0;
%print_Lagrangian_Beams(xLag,yLag,k_Beam,C,struct_name);


% Prints .target file!
k_Target = 5e6;
print_Lagrangian_Target_Pts(x1,k_Target,struct_name);

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
% FUNCTION: prints Vertex points to a file called rubberband.vertex
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Target_Pts(xLag,k_Target,struct_name)

    N = length(xLag);

    target_fid = fopen([struct_name '.target'], 'w');

    fprintf(target_fid, '%d\n', N );

    %Loops over all Lagrangian Pts.
    for s = 1:N
        if s > 284
            k_Target = 2e9;
        end
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
% FUNCTION: creates the Lagrangian structure geometry for PHASE 1
%           msg: "Hi KC!!"
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag] = give_Me_Immsersed_Boundary_Geometry_1(Lx,Nx,ds)

% The immsersed structure message #1: Hi KC!! %
len = Lx/5;

xC = -Lx/4;
yC = -Lx/2;
[xH,yH] = give_Me_The_Letter_Please(ds,len,'H',xC,yC);
[xi,yi] = give_Me_The_Letter_Please(ds,len,'i',xC-0.1,yC);
[xK,yK] = give_Me_The_Letter_Please(ds,len,'K',xC-0.3,yC);
[xCap,yCap] = give_Me_The_Letter_Please(ds,len,'C',xC-0.46,yC);
[xEx1,yEx1] = give_Me_The_Letter_Please(ds,1.2*len,'!',xC-0.54,yC);
[xEx2,yEx2] = give_Me_The_Letter_Please(ds,1.2*len,'!',xC-0.58,yC);

xBL2 =  0:ds:Lx/4+3*ds;  xBL2 = xBL2 + 3/8*Lx - 2*ds;
yTL2 =  0.75*ones(1,length(xBL2));
yBL2 =  0.25*ones(1,length(xBL2));

xBL = 0:ds:Lx/2+2*ds; xBL = xBL + Lx/4 - ds;
yTL = 0.82*ones(1,length(xBL));
yBL = 0.18*ones(1,length(xBL));

xBL3 = 0:ds:Lx/8+ds; xBL3 = xBL3 + 7/16*Lx -ds;
yBL3 = 0.9*ones(1,length(xBL3));
yTL3 = 0.1*ones(1,length(xBL3));


xLag = [xH xi xK xCap xEx1 xEx2 xBL xBL2 xBL3 xBL(1:end-1) xBL2 xBL3];
yLag = [yH yi yK yCap yEx1 yEx2 yBL yBL2 yBL3 yTL(1:end-1) yTL2 yTL3];

%plot(xH,yH,'r*'); hold on;
%plot(xi,yi,'r*'); hold on;
%plot(xK,yK,'r*'); hold on;
%plot(xCap,yCap,'r*'); hold on;
%plot(xEx1,yEx1,'r*'); hold on;
%plot(xEx2,yEx2,'r*'); hold on;
%axis([0 Lx 0 Lx]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates the Lagrangian structure geometry for PHASE 2
%           msg: "Would you like to ..."
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag] = give_Me_Immsersed_Boundary_Geometry_2(Lx,Nx,ds)

% The immsersed structure message #1: Hi KC!! %
len = Lx/8;

xC = -Lx/6;
yC = -2*Lx/3;
[xW,yW] = give_Me_The_Letter_Please(ds,len,'W',xC,yC);
[xo,yo] = give_Me_The_Letter_Please(ds,len,'o',xC-0.1,yC);
[xu,yu] = give_Me_The_Letter_Please(ds,len,'u',xC-0.18,yC);
[xl,yl] = give_Me_The_Letter_Please(ds,len,'l',xC-0.27,yC);
[xd,yd] = give_Me_The_Letter_Please(ds,1.2*len,'d',xC-0.29,yC);

[xy,yy] = give_Me_The_Letter_Please(ds,1.2*len,'y',xC-0.47,yC);
[xo2,yo2] = give_Me_The_Letter_Please(ds,len,'o',xC-0.55,yC);
[xu2,yu2] = give_Me_The_Letter_Please(ds,len,'u',xC-0.63,yC);

xC = -Lx/4;
yC = -1*Lx/4 - Lx/6;
[xl2,yl2] = give_Me_The_Letter_Please(ds,len,'l',xC,yC);
[xi,yi] = give_Me_The_Letter_Please(ds,len,'i',xC-0.01,yC);
[xk,yk] = give_Me_The_Letter_Please(ds,len,'k',xC-0.07,yC);
[xe,ye] = give_Me_The_Letter_Please(ds,len,'e',xC-0.13,yC);

[xt,yt] = give_Me_The_Letter_Please(ds,len,'t',xC-0.27,yC);
[xo3,yo3] = give_Me_The_Letter_Please(ds,len,'o',xC-0.32,yC);

yC = -1*Lx/5 - Lx/6;
[xo4,yo4] = give_Me_The_Letter_Please(ds/3,len/6,'o',xC-0.44,yC);
[xo5,yo5] = give_Me_The_Letter_Please(ds/3,len/6,'o',xC-0.48,yC);
[xo6,yo6] = give_Me_The_Letter_Please(ds/3,len/6,'o',xC-0.52,yC);

xBL = 0:ds:Lx/2.5; xBL = xBL + 3/10*Lx;
yTL = 0.84*ones(1,length(xBL));
yBL = 0.18*ones(1,length(xBL));


xLag = [xW xo xu xl xd xy xo2 xu2 xl2 xi xk xe xt xo3 xo4 xo5 xo6 xBL xBL];
yLag = [yW yo yu yl yd yy yo2 yu2 yl2 yi yk ye yt yo3 yo4 yo5 yo6 yTL yBL];

% plot(xW,yW,'r*'); hold on;
% plot(xo,yo,'r*'); hold on;
% plot(xu,yu,'r*'); hold on;
% plot(xl,yl,'r*'); hold on;
% plot(xd,yd,'r*'); hold on;
% 
% plot(xy,yy,'r*'); hold on;
% plot(xo2,yo2,'r*'); hold on;
% plot(xu2,yu2,'r*'); hold on;
% 
% plot(xl2,yl2,'r*'); hold on;
% plot(xi,yi,'r*'); hold on;
% plot(xk,yk,'r*'); hold on;
% plot(xe,ye,'r*'); hold on;
% 
% plot(xt,yt,'r*'); hold on;
% plot(xo3,yo3,'r*'); hold on;
% 
% plot(xo4,yo4,'r*'); hold on;
% plot(xo5,yo5,'r*'); hold on;
% plot(xo6,yo6,'r*'); hold on;
% axis([0 Lx 0 Lx]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates the Lagrangian structure geometry for PHASE 2
%           msg: "Would you like to ..."
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag] = give_Me_Immsersed_Boundary_Geometry_3(Lx,Nx,ds)

% The immsersed structure message #1: Hi KC!! %
len = Lx/9;

xC = -Lx/8;
yC = -3*Lx/4;
shift = 0.2;

[xo4,yo4] = give_Me_The_Letter_Please(ds/3,len/6,'o',xC,yC+0.09);
[xo5,yo5] = give_Me_The_Letter_Please(ds/3,len/6,'o',xC-0.04,yC+0.09);
[xo6,yo6] = give_Me_The_Letter_Please(ds/3,len/6,'o',xC-0.08,yC+0.09);

[xg,yg] = give_Me_The_Letter_Please(ds,len,'g',xC-shift,yC);
[xo,yo] = give_Me_The_Letter_Please(ds,len,'o',xC-0.1-shift,yC);

[xo2,yo2] = give_Me_The_Letter_Please(ds,1.2*len,'o',xC-0.25-shift,yC);
[xn,yn] = give_Me_The_Letter_Please(ds,1.2*len,'n',xC-0.34-shift,yC);

[xa,ya] = give_Me_The_Letter_Please(ds,len,'a',xC-0.48-shift,yC);


xC = -Lx/5;
yC = -Lx/2;
[xd,yd] = give_Me_The_Letter_Please(ds,len,'d',xC,yC);
[xa2,ya2] = give_Me_The_Letter_Please(ds,len,'a',xC-0.08,yC);
[xt,yt] = give_Me_The_Letter_Please(ds,len,'t',xC-0.16,yC);
[xe,ye] = give_Me_The_Letter_Please(ds,len,'e',xC-0.22,yC);

[xw,yw] = give_Me_The_Letter_Please(ds,1.1*len,'w',xC-0.37,yC);
[xi,yi] = give_Me_The_Letter_Please(ds,len,'i',xC-0.45,yC);
[xt2,yt2] = give_Me_The_Letter_Please(ds,len,'t',xC-0.505,yC);
[xh,yh] = give_Me_The_Letter_Please(ds,len,'h',xC-0.58,yC);

xC = -Lx/2.5;
yC = -Lx/4;
[xm,ym] = give_Me_The_Letter_Please(ds,1.1*len,'m',xC,yC);
[xe2,ye2] = give_Me_The_Letter_Please(ds,len,'e',xC-0.12,yC);
[xQU,yQU] = give_Me_The_Letter_Please(ds,1.2*len,'?',xC-0.22,yC);

xLag = [xo4 xo5 xo6 xg xo xo2 xn xa xd xa2 xt xe xw xi xt2 xh xm xe2 xQU];
yLag = [yo4 yo5 yo6 yg yo yo2 yn ya yd ya2 yt ye yw yi yt2 yh ym ye2 yQU];


% plot(xo4,yo4,'r*'); hold on;
% plot(xo5,yo5,'r*'); hold on;
% plot(xo6,yo6,'r*'); hold on;
% 
% plot(xg,yg,'r*'); hold on;
% plot(xo,yo,'r*'); hold on;
% 
% plot(xo2,yo2,'r*'); hold on;
% plot(xn,yn,'r*'); hold on;
% 
% plot(xa,ya,'r*'); hold on;
% 
% plot(xd,yd,'r*'); hold on;
% plot(xa2,ya2,'r*'); hold on;
% plot(xt,yt,'r*'); hold on;
% plot(xe,ye,'r*'); hold on;
% 
% plot(xw,yw,'r*'); hold on;
% plot(xi,yi,'r*'); hold on;
% plot(xt2,yt2,'r*'); hold on;
% plot(xh,yh,'r*'); hold on;
% 
% plot(xm,ym,'r*'); hold on;
% plot(xe2,ye2,'r*'); hold on;
% plot(xQU,yQU,'r*'); hold on;
% 
% axis([0 Lx 0 Lx]);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%s%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates the Lagrangian structure geometry for the Heart
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag] = give_Me_Immsersed_Boundary_Geometry_Heart(Lx,Nx,frac,Nb)

% The immsersed structure is a heart %
ds = 2*pi/(Nb-1);
tVec = 0:ds:2*pi;

for i=1:length(tVec)
    t = tVec(i);
    xLag(i)	=	16*sin(t)^3;
    yLag(i)	=	13*cos(t)-5*cos(2*t)-2*cos(3*t)-cos(4*t);
end

xLag = frac*xLag + 0.5;
yLag = frac*yLag + 0.5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: creates the Lagrangian structure geometry for PHASE 2 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag] = give_Me_Immsersed_Boundary_Geometry_4(Lx,Nx,frac,Nb)

% The immsersed structure is a heart %
ds = 2*pi/(Nb-1);
tVec = 0:ds:2*pi;
for i=1:length(tVec)
    t = tVec(i);
    r(i) = 1 - sin(t);
    xLag(i)	=	r(i)*cos(t);
    yLag(i)	=	r(i)*sin(t);
end

xLag = frac*xLag;
yLag = frac*yLag;

minY = min(yLag);
yLag = yLag - 3*minY/8;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: take out Lagrangian Pts.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x1,y1,x2,y2] = please_Take_Out_Points(x1,y1,x2,y2)

n1 = 50;
n2 = 70;

x1 = [x1(1:n1) x1(n2:end)];
y1 = [y1(1:n1) y1(n2:end)];

x2 = [x2(1:n1) x2(n2:end)];
y2 = [y2(1:n1) y2(n2:end)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints all Vertices to File
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function please_Print_Vertices_To_File(X1,Y1,X2,Y2,X3,Y3)

fileID = fopen('All_Positions.txt','w');
for j=1:length(X1)
    fprintf(fileID,'%1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e\n', X1(j),Y1(j),X2(j),Y2(j),X3(j),Y3(j),X3(j),Y3(j) );
end
fclose(fileID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: give me the letter!!!!!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,y] = give_Me_The_Letter_Please(ds,len,letter,xC,yC)

if strcmp(letter,'a')
    r = len/4.1;
    dt = ds/r;
    theta = 0:dt:2*pi;
    for i=1:length(theta)
       xA(i) = r*cos(theta(i)) - xC;
       yA(i) = r*sin(theta(i)) - yC - len/4;
    end
    
    yLine = [-r:ds:r r]; yLine = yLine - yC - len/4;
    xLine = r*ones(1,length(yLine)) - (xC-ds/6);
    
    x = [xA xLine];
    y = [yA yLine];

elseif strcmp(letter,'C')
    
    r = len/2;
        dt = ds/r;
    theta = pi/2:dt:pi-0.25*dt;
    theta = [theta pi -theta];
    for i=1:length(theta)
       x(i) = 3*r/4*cos(theta(i)) - xC;
       y(i) = r*sin(theta(i)) - yC;
    end
    
    xTB = 0:3*ds/4:len/6;
    xTB = xTB - xC;
    yT = r*ones(1,length(xTB))-yC;
    yB = -r*ones(1,length(xTB))-yC;
    x = [xTB x xTB];
    y = [yT y yB];
    
elseif strcmp(letter,'d')
    
    r = len/4.1;
    dt = ds/r;
    theta = 0:dt:2*pi;
    for i=1:length(theta)
       xA(i) = r*cos(theta(i)) - xC;
       yA(i) = r*sin(theta(i)) - yC - len/4;
    end
    
    yLine = [-len/2:ds:len/2.25 r]; yLine = yLine - yC;
    xLine = r*ones(1,length(yLine)) - (xC-ds/6);
    
    x = [xA xLine];
    y = [yA yLine];
    
 elseif strcmp(letter,'e')
        
    r = len/4.1;
    yV = [-r:ds:r r]; 
    xH = yV - xC;
    yH = zeros(1,length(xH)) - yC;
    yH2 = yH - len/4;
    yH3 = yH - len/2;
    
    yV = yV - yC - len/4;
    xV = r*ones(1,length(yV));
    
    yV3 = 0:-ds:-len/4;
    yV3 = yV3 - yC;
    xV3 = r*ones(1,length(yV3));
        
    x = [-xV-(xC-ds/6) xH xV3-(xC-ds/6) xH xH];
    y = [yV yH yV3 yH2 yH3];
           
    
elseif strcmp(letter,'g')
    
    r = len/4.1;
    dt = ds/r;
    theta = 0:dt:2*pi;
    for i=1:length(theta)
       xA(i) = r*cos(theta(i)) - xC;
       yA(i) = r*sin(theta(i)) - yC - len/4;
    end
    
    yLine = [-len/2:ds:len/2.25 r]; yLine = yLine - yC - len/3;
    xLine = r*ones(1,length(yLine)) - (xC-ds/6);
    
    yV = [-r:ds:r r]; 
    xH = yV - xC;
    yH = zeros(1,length(xH)) - yC - 0.85*len;
    
    x = [xA xLine xH];
    y = [yA yLine yH];
    
elseif strcmp(letter,'h')

   yLine = -len/2:ds:0-ds/4;
   yLine = [yLine 0 -yLine];
   yLine1 = yLine - yC;
   xLine1 = 1/2*len/2*ones(1,length(yLine1));
   
   yLine2 = -len/2:ds:0-ds/4;
   yLine2 = [yLine2 0];
   yLine2 = yLine2 - yC;
   xLine2 = 1/2*len/2*ones(1,length(yLine2));
   
   xH = (-1/2*len/2+ds/2:ds:1/2*len/2-ds/4);
   xH = xH - xC;
   yH = zeros(1,length(xH)) - yC;
   
   x = [-xLine1-xC xH xLine2-xC];
   y = [yLine1 yH yLine2];    
   
elseif strcmp(letter,'H')
   yLine = -len/2:ds:0-ds/4;
   yLine = [yLine 0 -yLine];
   yLine = yLine - yC;
   xLine = 1/2*len/2*ones(1,length(yLine));
   xH = (-1/2*len/2+ds/2:ds:1/2*len/2-ds/4);
   xH = xH - xC;
   yH = zeros(1,length(xH)) - yC;
   
   x = [-xLine-xC xH xLine-xC];
   y = [yLine yH yLine];
   
elseif strcmp(letter,'i')
    
   yLine = -len/2:ds:0-ds/4;
   yLine = yLine - yC;
   xLine = zeros(1,length(yLine)) - xC; 
   
    r = len/21;
    dt = ds/(2*r);
    theta = 0:dt:2*pi;
    for i=1:length(theta)
       xDot(i) = r*cos(theta(i)) - xC;
       yDot(i) = r*sin(theta(i)) - yC + len/7;
    end
    
    x = [xLine xDot];
    y = [yLine yDot];
    

elseif strcmp(letter,'k')
    
   yLine = -len/2:ds:len/2;
   yLine = yLine - yC;
   xLine = zeros(1,length(yLine)) -xC-len/3.9;

   XD = 0:ds:len/sqrt(2);
   YD = zeros(1,length(XD));
   
   %rotate!
   for i=1:length(XD)
      xDt(i) = XD(i)*cos(pi/4)-YD(i)*sin(pi/4);
      yDt(i) = XD(i)*sin(pi/4)+YD(i)*cos(pi/4);
      
      xDb(i) = XD(i)*cos(pi/4)-YD(i)*sin(-pi/4);
      yDb(i) = XD(i)*sin(-pi/4)+YD(i)*cos(pi/4);
   end
   xDt = xDt - xC - len/4;
   yDt = yDt - yC;
   
   xDb = xDb - xC - len/4;
   yDb = yDb - yC;
   
   x = [xLine(1:ceil(0.75*end)) xDt(1:ceil(2*end/3)) xDb(1:ceil(2*end/3))];
   y = [yLine(1:ceil(0.75*end)) yDt(1:ceil(2*end/3))-len/8 yDb(1:ceil(2*end/3))-len/5.5];
    
elseif strcmp(letter,'K')
    
   yLine = -len/2:ds:len/2;
   yLine = yLine - yC;
   xLine = zeros(1,length(yLine)) -xC-len/3.9;

   XD = 0:ds:len/sqrt(2);
   YD = zeros(1,length(XD));
   
   %rotate!
   for i=1:length(XD)
      xDt(i) = XD(i)*cos(pi/4)-YD(i)*sin(pi/4);
      yDt(i) = XD(i)*sin(pi/4)+YD(i)*cos(pi/4);
      
      xDb(i) = XD(i)*cos(pi/4)-YD(i)*sin(-pi/4);
      yDb(i) = XD(i)*sin(-pi/4)+YD(i)*cos(pi/4);
   end
   xDt = xDt - xC - len/4;
   yDt = yDt - yC;
   
   xDb = xDb - xC - len/4;
   yDb = yDb - yC;
   
   x = [xLine xDt xDb];
   y = [yLine yDt yDb];

 
elseif strcmp(letter,'l')
    
   yLine = -len/2:ds:len/2.25;
   yLine = yLine - yC;
   xLine = 1/2*len/2*ones(1,length(yLine));

   x = -xLine-xC;
   y = yLine;

elseif strcmp(letter,'m')
        
    r = len/2;
    xH = -r:ds:r;  xH = xH - xC;
    yH = zeros(1,length(xH)) - yC;
    yV = 0:-ds:-len/2; yV = yV - yC;
    xV = -r*ones(1,length(yV)) - xC;
    xV2 = zeros(1,length(yV)) - xC;
    xV3 = r*ones(1,length(yV)) - xC;
    
    yVS = 0:ds:len/12;
    yVS = yVS - yC;
    xVS = -r*ones(1,length(yVS)) - xC;
        
    x = [xVS xH xV xV2 xV3];
    y = [yVS yH yV yV yV];
    
elseif strcmp(letter,'n')
        
    r = len/4.1;
    yV = [-r:ds:r r]; 
    xH = yV - xC;
    yH = zeros(1,length(xH)) - yC;
    
    yV = yV - yC - len/4;
    xV = r*ones(1,length(yV));
    
    yV2 = 0:ds:len/12;
    yV2 = yV2 - yC;
    xV2 = r*ones(1,length(yV2));
        
    x = [-xV-(xC-ds/6) -xV2-(xC-ds/6) xH xV-(xC-ds/6)];
    y = [yV yV2 yH yV];
    
   
elseif strcmp(letter,'o')
       
    r = len/4.1;
    dt = ds/r;
    theta = 0:dt:2*pi;
    for i=1:length(theta)
       x(i) = r*cos(theta(i)) - xC;
       y(i) = r*sin(theta(i)) - yC - len/4;
    end
    
 elseif strcmp(letter,'r')
        
    r = len/4.1;
    yV = [-r:ds:r r]; 
    xH = yV - xC;
    yH = zeros(1,length(xH)) - yC;
    
    yV = yV - yC - len/4;
    xV = r*ones(1,length(yV));
    
    yV2 = 0:ds:len/12;
    yV2 = yV2 - yC;
    xV2 = r*ones(1,length(yV2));
    
    yV3 = 0:-ds:-len/10;
    yV3 = yV3 - yC;
    xV3 = r*ones(1,length(yV3));
        
    x = [-xV-(xC-ds/6) -xV2-(xC-ds/6) xH xV3-(xC-ds/6)];
    y = [yV yV2 yH yV3];
       
    
elseif strcmp(letter,'t')
    
   yLine = -len/2:ds:len/2.25;
   yLine = yLine - yC;
   xLine = zeros(1,length(yLine));

   xH = -len/4:ds:len/4;
   xH = xH-xC;
   yH = zeros(1,length(xH)) - yC + len/9;
   
   x = [xLine-xC xH];
   y = [yLine yH];
   
    
elseif strcmp(letter,'u')
        
    r = len/4.1;
    yV = [-r:ds:r r]; 
    xH = yV - xC;
    yH = zeros(1,length(xH)) - yC - len/2;
    
    yV = yV - yC - len/4;
    xV = r*ones(1,length(yV));
    
    x = [-xV-(xC-ds/6) xH xV-(xC-ds/6)];
    y = [yV yH yV];

elseif strcmp(letter,'v')
    
   XD = 0:ds:sqrt(5)*len/4;
   YD = zeros(1,length(XD));
   
   %rotate!
   ang = atan(2);
   for i=1:length(XD)
      xDr(i) = XD(i)*cos(ang)-YD(i)*sin(ang);
      yDr(i) = XD(i)*sin(ang)+YD(i)*cos(ang);
      
      xDL(i) = XD(i)*cos(pi-ang)-YD(i)*sin(pi-ang);
      yDL(i) = XD(i)*sin(pi-ang)+YD(i)*cos(pi-ang);
   end
   xDr = xDr - xC;
   yDr = yDr - yC - len/2;
 
   xDL = xDL - xC;
   yDL = yDL - yC - len/2;
 
   x = [xDL xDr(2:end)];
   y = [yDL yDr(2:end)];
   
elseif strcmp(letter,'w')
    
   XD = 0:ds:sqrt(5)*len/4;
   YD = zeros(1,length(XD));
   
   %rotate!
   ang = atan(2);
   for i=1:length(XD)
      xDr(i) = XD(i)*cos(ang)-YD(i)*sin(ang);
      yDr(i) = XD(i)*sin(ang)+YD(i)*cos(ang);
      
      xDL(i) = XD(i)*cos(pi-ang)-YD(i)*sin(pi-ang);
      yDL(i) = XD(i)*sin(pi-ang)+YD(i)*cos(pi-ang);
   end
   xDr = xDr - xC - len/4;
   yDr = yDr - yC - len/2;
 
   xDL = xDL - xC - len/4;
   yDL = yDL - yC - len/2;
 
   x = [xDL xDr(2:end) xDL+len/2 xDr(2:end)+len/2];
   y = [yDL yDr(2:end) yDL yDr(2:end)];   
   
elseif strcmp(letter,'W')
    
   XD = 0:ds:sqrt(5)*len/2;
   YD = zeros(1,length(XD));
   
   %rotate!
   ang = atan(2);
   for i=1:length(XD)
      xDr(i) = XD(i)*cos(ang)-YD(i)*sin(ang);
      yDr(i) = XD(i)*sin(ang)+YD(i)*cos(ang);
      
      xDL(i) = XD(i)*cos(pi-ang)-YD(i)*sin(pi-ang);
      yDL(i) = XD(i)*sin(pi-ang)+YD(i)*cos(pi-ang);
   end
   xDr = xDr - xC - len/4;
   yDr = yDr - yC - len/2;
 
   xDL = xDL - xC - len/4;
   yDL = yDL - yC - len/2;
 
   x = [xDL xDr(2:floor(end/2)) xDL(1:floor(end/2))+len/2 xDr(2:end)+len/2];
   y = [yDL yDr(2:floor(end/2)) yDL(1:floor(end/2))       yDr(2:end)];     
   
elseif strcmp(letter,'y')
    
   XD = 0:ds:sqrt(5)*len/4;
   YD = zeros(1,length(XD));
   
   %rotate!
   ang = atan(2);
   for i=1:length(XD)
      xDr(i) = XD(i)*cos(ang)-YD(i)*sin(ang);
      yDr(i) = XD(i)*sin(ang)+YD(i)*cos(ang);
      
      xDL(i) = XD(i)*cos(pi-ang)-YD(i)*sin(pi-ang);
      yDL(i) = XD(i)*sin(pi-ang)+YD(i)*cos(pi-ang);
      
      xB(i) = XD(i)*cos(pi+ang)-YD(i)*sin(pi+ang);
      yB(i) = XD(i)*sin(pi+ang)+YD(i)*cos(pi+ang);
   end
   xDr = xDr - xC;
   yDr = yDr - yC - len/2;
 
   xDL = xDL - xC;
   yDL = yDL - yC - len/2;
   
   xB = xB - xC;
   yB = yB - yC - len/2;
 
   x = [xDL xDr(2:end) xB(2:end)];
   y = [yDL yDr(2:end) yB(2:end)];
      
   
elseif strcmp(letter,'!')
    
   yLine = -len/4:ds:len/2;
   yLine = yLine - yC;
   xLine = zeros(1,length(yLine)) - xC; 
   
    r = len/21;
    dt = ds/(2*r);
    theta = 0:dt:2*pi;
    for i=1:length(theta)
       xDot(i) = r*cos(theta(i)) - xC;
       yDot(i) = r*sin(theta(i)) - yC - len/2.1;
    end
    
    x = [xLine xDot];
    y = [yLine yDot];  
    
elseif strcmp(letter,'?')
    
    %xT = 0:ds:len/8-ds; xT = xT - xC - len/4;
    %yT = zeros(1,length(xT)) - yC + len/2;

    r = len/4;
    dt = ds/r;
    theta = 1.25*pi/2:-dt:-pi/2;
    for i=1:length(theta)
       xL(i) = 1.5*r*cos(theta(i)) - xC - len/8;
       yL(i) = r*sin(theta(i)) - yC + len/4;
    end
    
    yV = -ds:-ds:-len/3.25;
    xV = zeros(1,length(yV)) - xC - len/8;
    
    r = len/21;
    dt = ds/(2*r);
    theta = 0:dt:2*pi;
    for i=1:length(theta)
       xDot(i) = r*cos(theta(i)) - xC -len/8;
       yDot(i) = r*sin(theta(i)) - yC - len/2.1;
    end
    
    x = [xL xV xDot];
    y = [yL yV yDot];  
 
end



