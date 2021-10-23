function give_Me_Jelly()

L = 3;           % Length of Computational Domain
ds = L/(2*600);  % Makes Lagrangian Spaced based on 1024x1024 finest grid

frac = 0.5;

% Makes Geometry for EQUILIBRIUM State %
%
a1 = 0.63;  %Semi-major axis length for position ONE (Expanded State)
b1 = 0.45;  %Semi-minor axis length for position ONE (Expanded State)
a1 = frac*a1; b1 = frac*b1;
a1
b1
[X1 Y1 Ninfo arcL] = make_Geometry_Position_ONE(ds,a1,b1);

% Makes Geometry for OVER-CONTRACTED State for RIGHT side %
% 
a2 = 0.325; %Semi-major axis length for position TWO (Contracted State)
b2 = 0.6;   %GUESS: Semi-minor axis length for position TWO (Contracted State)
a2 = frac*a2; b2 = frac*b2;
[X2 Y2 b2] = make_Geometry_Position_TWO(a2,b2,b1,Ninfo,arcL);

% Makes Geometry for LESS-CONTRACTED State for LEFT side %
%
a3 = 0.4;   %Semi-major axis length for position TWO (Contracted State)
b3 = 0.7;   %GUESS: Semi-minor axis length for position TWO (Contracted State)
a3 = frac*a3; b3 = frac*b3;
a3
b3
[X3 Y3 b3] = make_Geometry_Position_TWO(a3,b3,b1,Ninfo,arcL);


% Makes Geometry for OVER-EXPANDED State for LEFT SIDE %
%
a4 = 0.75;  %Semi-major axis length for position TWO (Contracted State)
b4 = 0.3;   %GUESS: Semi-minor axis length for position TWO (Contracted State)
a4 = frac*a4; b4 = frac*b4;
[X4 Y4 b4] = make_Geometry_Position_TWO(a4,b4,b1,Ninfo,arcL);

% Translate Points from [-L/2,L/2]x[-L/2,L/2] -> [0,L]x[0,L]
X1 = X1 + L/2; Y1 = Y1 + L/2; 
X2 = X2 + L/2; Y2 = Y2 + L/2; 
X3 = X3 + L/2; Y3 = Y3 + L/2; 
X4 = X4 + L/2; Y4 = Y4 + L/2; 

nArm = (length(X1)-1)/2;

% Plots Geometry for Both Phases %
val = frac*0.75;
plot(X1,Y1,'*',X2(1:nArm),Y2(1:nArm),'r*',X3(nArm+2:end),Y3(nArm+2:end),'g*',X4(nArm+2:end),Y4(nArm+2:end),'m*'); hold on;
axis([0 L 0 L]);
legend('1','2','3','4');

figure(2)
plot(X2,Y2,'r*'); hold on;
plot(X1,Y1,'b*'); hold on;

% Prints Useful Info about Geometry %
fprintf('\n\n');
fprintf('Number of Total Geometry Pts: %d\n',Ninfo(2));
fprintf('Number of Geometry Pts in EACH Arm: %d\n',Ninfo(1));
fprintf('Arc length of bell (one side): %d\n',arcL);
fprintf('\n');
fprintf('(X1,Y1): Equilibrium State\n');
fprintf('(X2,Y2): OVER-Contracted State for RIGHT Side\n');
fprintf('(X3,Y3): Contracted State for LEFT Side (slightly less than right side)\n');
fprintf('(X4,Y4): OVER-EXPANDED State for LEFT Side\n');
fprintf('\n');
fprintf('Idea: \n');
fprintf('RHS contracts further than LHS, rests & returns to equil. state after LHS begins expanding\n');
fprintf('LHS contracts less than RHS, then expands further than equil. state and then returns to equil. state\n');
fprintf('\n');
fprintf('PUT FOLLOWING IN UPDATE_TARGET_POINTS.C:\n');
fprintf('numPts (# of Total Geometry Pts): %d\n',Ninfo(2));
fprintf('NL_Start (FIRST Point in LEFT Arm): %d\n\n\n',Ninfo(1)+1);


% Prints INPUT files %
print_Vertex_Pts(X1,Y1,Ninfo);
target_force = 5e6;
print_Target_Pts(target_force,Ninfo);


% Prints all vectors to .TXT files %
%print_Text_Files(Ninfo,X1,Y1,X2,Y2,X3,Y3,X4,Y4)

% Prints 2-Position Right/Left Side vectors to .TXT files %
print_2_Position_Text_File(Ninfo,X1,Y1,X2,Y2,X3,Y3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: Prints 2 Position Phases:
%           XY(:,1) = original position x-Values
%           XY(:,3) = original position y-Values
%           XY(2:N+1,2) = 2nd position x-Values for RIGHT bell
%           XY(N+2:end,2) = 2nd position x-Values for LEFT bell
%           XY(2:N+1,4) = 2nd position y-Values for RIGHT bell
%           XY(N+2:end,4) = 2nd position y-Values for LEFT bell
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_2_Position_Text_File(Ninfo,X1,Y1,X2,Y2,X3,Y3)


N = Ninfo(1);    % Number of pts. in ONE arm
NTot = Ninfo(2); % Total # of Lag Pts. on Jelly Bell

% Initialization
XY = zeros(NTot,4); 

% Construct Matrices For Storing Points
XY(:,1) = X1(1:NTot); XY(:,2) = X1(1:NTot); 
XY(:,3) = Y1(1:NTot); XY(:,4) = Y1(1:NTot); 

% 2nd Phase Right Arm
XY(2:N+1,2) = X2(2:N+1);
XY(2:N+1,4) = Y2(2:N+1);


% 2nd Phase Left Arm
XY(N+2:end,2) = X2(N+2:end);
XY(N+2:end,4) = Y2(N+2:end);

% PLOT TEST! %
%figure(3)
%plot(XY(:,1),XY(:,3),'k*'); hold on;
%plot(XY(:,2),XY(:,4),'g*'); hold on;


% PRINT x/y-Information %

fileID_XY = fopen('XY_2Pos.txt','w');

for i=1:length( XY(:,1) )
    fprintf(fileID_XY,'%1.16e %1.16e %1.16e %1.16e\n',XY(i,1),XY(i,2),XY(i,3),XY(i,4));
end

fclose(fileID_XY);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to PRINT all TEXT files for both PHASES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Text_Files(Ninfo,X1,Y1,X2,Y2,X3,Y3,X4,Y4)

N = Ninfo(1); %Number of pts. in ONE arm

% x-Information %

fileID_x1 = fopen('x1.txt','w');
fprintf(fileID_x1,'%1.16e\n',X1);
fclose(fileID_x1);

fileID_x2 = fopen('x2.txt','w');
fprintf(fileID_x2,'%1.16e\n',X2);
fclose(fileID_x2);

fileID_x3 = fopen('x3.txt','w');
fprintf(fileID_x3,'%1.16e\n',X3);
fclose(fileID_x3);

fileID_x4 = fopen('x4.txt','w');
fprintf(fileID_x4,'%1.16e\n',X4);
fclose(fileID_x4);

% y-Information %

fileID_y1 = fopen('y1.txt','w');
fprintf(fileID_y1,'%1.16e\n',Y1);
fclose(fileID_y1);

fileID_y2 = fopen('y2.txt','w');
fprintf(fileID_y2,'%1.16e\n',Y2);
fclose(fileID_y2);

fileID_y3 = fopen('y3.txt','w');
fprintf(fileID_y3,'%1.16e\n',Y3);
fclose(fileID_y3);

fileID_y4 = fopen('y4.txt','w');
fprintf(fileID_y4,'%1.16e\n',Y4);
fclose(fileID_y4);


% xR1 = X1(2:2+N-1); % x-Pts. on RHS of center point for Position ONE
% xL1 = X1(2+N:end); % x-Pts. on LHS of center point for Position ONE
% 
% xR2 = X2(2:2+N-1); % x-Pts. on RHS of center point for Position TWO
% xL2 = X2(2+N:end); % x-Pts. on LHS of center point for Position TWO
% 
% yR1 = Y1(2:2+N-1); % y-Pts. on RHS of center point for Position ONE
% yL1 = Y1(2+N:end); % y-Pts. on RHS of center point for Position ONE
% 
% yR2 = Y2(2:2+N-1); % y-Pts. on LHS of center point for Position TWO
% yL2 = Y2(2+N:end); % y-Pts. on LHS of center point for Position TWO
% 
% 
% fileID_xR1 = fopen('xR_1.txt','w');
% fprintf(fileID_xR1,'%1.16e\n',xR1);
% fclose(fileID_xR1);
% 
% fileID_xR2 = fopen('xR_2.txt','w');
% fprintf(fileID_xR2,'%1.16e\n',xR2);
% fclose(fileID_xR2);
% 
% fileID_xL1 = fopen('xL_1.txt','w');
% fprintf(fileID_xL1,'%1.16e\n',xL1);
% fclose(fileID_xL1);
% 
% fileID_xL2 = fopen('xL_2.txt','w');
% fprintf(fileID_xL2,'%1.16e\n',xL2);
% fclose(fileID_xL2);
% 
% 
% fileID_yR1 = fopen('yR_1.txt','w');
% fprintf(fileID_yR1,'%1.16e\n',yR1);
% fclose(fileID_yR1);
% 
% fileID_yR2 = fopen('yR_2.txt','w');
% fprintf(fileID_yR2,'%1.16e\n',yR2);
% fclose(fileID_yR2);
% 
% fileID_yL1 = fopen('yL_1.txt','w');
% fprintf(fileID_yL1,'%1.16e\n',yL1);
% fclose(fileID_yL1);
% 
% fileID_yL2 = fopen('yL_2.txt','w');
% fprintf(fileID_yL2,'%1.16e\n',yL2);
% fclose(fileID_yL2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Makes the Oblate Geometry (via ellipses)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X1 Y1 Ninfo arcL] = make_Geometry_Position_ONE(ds,a,b)

%a: Semi-major axis (along horizontal)
%b: Semi-minor axis (along vertical)

%Finding elliptical arc length
h = ( (a-b)/(a+b) )^2;
arcL = (1/4)* pi*(a+b)*(1 + 3*h/ (10 + sqrt(4-3*h) ) ); % (1/4)th of the entire ELLIPTICAL ARC

n = 1;        %counting for array
X(n) = 0;     %first x-value at (0,-b)
Y(n) = -b;    %first y-value at (0,-b)
tprev = -pi/2;  %initial angle 

tol = ds / 10; %error tolerance for root-finding algorithm

maxy = max(a,b);

    while tprev < 0 
    
        n = n+1;
        tnext = tprev + ds/maxy; %Guess for next angle value
        
        if tprev < -pi/4
            tfar =  0.0;%abs(2*tprev);      %Far guess
        else
            tfar = pi/4;
        end
            
        %initiating guess for bisection-algorithm
        xn = a*cos(tnext);
        yn = b*sin(tnext);
        errSign = ( ds - sqrt( (xn-X(n-1))^2 + (yn-Y(n-1))^2 ) );
        err = abs(errSign);
    
        %Bisection algorithm to make points equally spaced
        while ( err > tol )
        
            if errSign < 0
                tfar = tnext;
                tnext = (tnext+tprev)/2;
            elseif errSign > 0
                tprev = tnext;
                tnext = (tnext+tfar)/2;
            end
            
            if tnext<tprev
                fprintf('NOT CONVERGING AT ANGLE %d\n',tprev)
                break;  
            end
            
            xn = a*cos(tnext);
            yn = b*sin(tnext);
            errSign = ( ds - sqrt( (xn-X(n-1))^2 + (yn-Y(n-1))^2 ) );
            err = abs(errSign);
        end
        
        X(n) = xn;   %Store X-value
        Y(n) = yn;   %Store Y-value
        tprev = tnext; %Update previous angle
    
        %plot(xn,yn,'*'); hold on;

        
    end
    
 N = length(X);   
 X(N+1) = a; %Add one more pt.
 Y(N+1) = 0; %Add one more pt.
  
 X1 = [X -X(2:end)]; %Reflects geometry 
 Y1 = -[Y  Y(2:end)]; %Y-values for reflected X-value points
 
 Ninfo(1) = N;          % # of pts. along one arm (NOT counting center pt.)
 Ninfo(2) = length(X1); % Total # of pts. in geometry

 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Makes the Oblate Geometry (via ellipses)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X2 Y2 b] = make_Geometry_Position_TWO(a,bguess,b1,Ninfo,L)


%Find correct b to give same length arm
b = Newtons_Give_Me_Semi_Minor(a,bguess,L);


%Finding elliptical arc length
h = ( (a-b)/(a+b) )^2;
arcL = (1/4)* pi*(a+b)*(1 + 3*h/ (10 + sqrt(4-3*h) ) ); % (1/4)th of the entire ELLIPTICAL ARC

ds2 = arcL / Ninfo(1); %Ensures same number of pts on arm for position ONE and TWO

[X Y Ninfo] = make_Elliptic_Geometry_Position_TWO(ds2,a,b,Ninfo);

dV = b1-b;   %Makes sure first pt is the same as in position ONE

Y2 = -(Y-dV);%Minus sign to put jelly into right orientation
X2 = X;

%Ninfo
%plot(X2,Y2,'r*'); hold on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Makes the Oblate Geometry (via ellipses)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function b = Newtons_Give_Me_Semi_Minor(a,bn,L)

tol = 1e-6;

bprev = 0.2*bn; %"Previous Guess" (left side of guess)
bnext = bn;     % First "Guess" at b
bfar =  1.5*bn; %"Far Guess" (right side of guess)

%Finding elliptical arc length
h = ( (a-bn)/(a+bn) )^2;
arcL = (1/4)* pi*(a+bn)*(1 + 3*h/ (10 + sqrt(4-3*h) ) ); % (1/4)th of the entire ELLIPTICAL ARC

errSign = L - arcL;
err = abs( errSign );


%Bisection algorithm to make points equally spaced
while ( err > tol )
        
    if errSign < 0
        bfar = bnext;
        bnext = (bnext+bprev)/2;
    elseif errSign > 0
        bprev = bnext;
        bnext = (bnext+bfar)/2;
    end
            
    %if bnext<bprev
    %    fprintf('NOT CONVERGING AT ANGLE %d\n',tprev)
    %    break;  
    %end
         
    %Finding elliptical arc length
    h = ( (a-bnext)/(a+bnext) )^2;
    arcL = (1/4)* pi*(a+bnext)*(1 + 3*h/ (10 + sqrt(4-3*h) ) ); % (1/4)th of the entire ELLIPTICAL ARC
    
    errSign = ( L - arcL );
    err = abs(errSign);
end

b = bnext; %Assign value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Makes the Oblate Geometry (via ellipses)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X1 Y1 Ninfo] = make_Elliptic_Geometry_Position_TWO(ds,a,b,Ninfo)

%a: Semi-major axis (along horizontal)
%b: Semi-minor axis (along vertical)

n = 1;        %counting for array
X(n) = 0;     %first x-value at (0,-b)
Y(n) = -b;    %first y-value at (0,-b)
tprev = -pi/2;  %initial angle 

tol = ds / 10; %error tolerance for root-finding algorithm

maxy = max(a,b);

    while n < Ninfo(1)+1 
    
        n = n+1;
        tnext = tprev + ds/maxy; %Guess for next angle value
        
        if tprev < -pi/4
            tfar =  0.0;%abs(2*tprev);      %Far guess
        else
            tfar = pi/4;
        end
            
        %initiating guess for bisection-algorithm
        xn = a*cos(tnext);
        yn = b*sin(tnext);
        errSign = ( ds - sqrt( (xn-X(n-1))^2 + (yn-Y(n-1))^2 ) );
        err = abs(errSign);
    
        %Bisection algorithm to make points equally spaced
        while ( err > tol )
        
            if errSign < 0
                tfar = tnext;
                tnext = (tnext+tprev)/2;
            elseif errSign > 0
                tprev = tnext;
                tnext = (tnext+tfar)/2;
            end
            
            if tnext<tprev
                fprintf('NOT CONVERGING AT ANGLE %d\n',tprev)
                break;  
            end
            
            xn = a*cos(tnext);
            yn = b*sin(tnext);
            errSign = ( ds - sqrt( (xn-X(n-1))^2 + (yn-Y(n-1))^2 ) );
            err = abs(errSign);
        end
        
        X(n) = xn;   %Store X-value
        Y(n) = yn;   %Store Y-value
        tprev = tnext; %Update previous angle
    
        %plot(xn,yn,'*'); hold on;

        
    end
    
 N = length(X)-1;   
 %X(N+1) = a; %Add one more pt.
 %Y(N+1) = 0; %Add one more pt.
 
 %plot(X,Y,'*'); hold on;
 
 X1 = [X -X(2:end)]; %Reflects geometry 
 Y1 = [Y  Y(2:end)]; %Y-values for reflected X-value points
 
 Ninfo(1) = N;          % # of pts. along one arm (NOT counting center pt.)
 Ninfo(2) = length(X1); % Total # of pts. in geometry














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to print VERTEX points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Vertex_Pts(X,Y,Ninfo)

%Ninfo(1): # of pts. on one arm on Jellyfish
%Ninfo(2): # of total pts. in Jellyfish Geometry

Ntot = Ninfo(2);

vertex_fid = fopen('jellyfish.vertex', 'w');

% first line is the number of vertices in the file
fprintf(vertex_fid, '%d\n', Ntot);

for i=1:Ntot
    fprintf(vertex_fid, '%1.16e %1.16e\n', X(i), Y(i) );    
end
    
fclose(vertex_fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to print TARGET pts.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Target_Pts(target_force,Ninfo)

%Ninfo(1): # of pts. on one arm on Jellyfish
%Ninfo(2): # of total pts. in Jellyfish Geometry

Ntot = Ninfo(2);

% Write out the target point information
target_fid = fopen('jellyfish.target', 'w');

% First line is the number of target pts. in the file
fprintf(target_fid, '%d\n', Ntot);

for s=1:Ntot
        fprintf(target_fid, '%d %1.16e\n', s , target_force);
end

fclose(target_fid);