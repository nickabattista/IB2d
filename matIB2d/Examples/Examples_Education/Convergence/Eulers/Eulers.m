%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Nicholas A. Battista
% Institution: The College of New Jersey (TCNJ)
% Created: August 12, 2014
% Date of Last Revision: February 3, 2019
%
% PURPOSE: -To demonstrate a convergence study using the Euler Method
% 
% Specifics: 1. This function solves the following ODE:
%               dy/dt = f(t,y)
%               y(0) = y0
%               using Eulers Method. 
%
%            2. It then does a convergence study for various h values. 
%
% Note: It performs the convergence study for the known ODEs,
%
% dy/dt = y, with y(0) = 1; w/ exact solution is: y(t) = e^t
% and
% dy/dt = 2*pi*cos(2*pi*t), w/ y(0)=1 w/ exact solution y(t) = sin(2*pi*t)+1
%
% Inputs:   
%
% Returns:  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Eulers()

% Prints information about script
print_info();

y0 = 1;   %Initial Condition, y(0)=y0

tS = 0;   %Starting Time for Simulation
tE = 2.0; %End Time for Simulation

%Makes vector of time iterates for convergence study
NVec = give_Me_NumberOfGridPts(); 

%Allocate memory to storage vectors
hVec = zeros(1,length(NVec));
err_h = hVec;
Y_SOL = zeros(1e5,5);
ct = 1;

%For loop for convergence study
for j=1:length(NVec)
    
    %Begin counting time for integration
    tic
    
    N = NVec(j);    %# of time-steps
    h = (tE-tS)/N;  %time-step
    hVec(j) = h;    %stores time-step
    
    time=tS:h:tE;   
    yFOR = zeros(1,length(time));
    yFOR2= yFOR;
    
    %Performs the Euler Method Time-Stepping Scheme
    for i=1:length(time)
    
        if i==1
            yFOR(i) = y0;
            yFOR2(i)= y0;
        else
            yFOR(i) = yFOR(i-1) + h*f(time(i-1),yFOR(i-1),1);
            yFOR2(i)= yFOR2(i-1)+ h*f(time(i-1),yFOR(i-1),2);
        end
        
    end

    %Gives Error at each time-step
    err = compute_Error(time,yFOR,y0,1);
    err2= compute_Error(time,yFOR2,y0,2);

    %Computes Absolute Error ("Inf-Norm of Error")
    err_h(j)  = max( abs(err) );
    err_h2(j) = max( abs(err2) );

    %Stores time for computation
    timeV(j) = toc/2;
    
    %Stores numerical solution / info for plotting soln.
    if ( ( (mod(j,3)==0) && (j<10) ) || (j==11) || (j==14) )
        Y_SOL(1:N+1,ct) = yFOR;
        Y_SOL2(1:N+1,ct)= yFOR2;
        hVecPlot(ct) = h;
        NVecPlot(ct) = N;
        errMat(1:N+1,ct) = err;
        errMat2(1:N+1,ct) = err2;
        ct = ct+1;
    end
    clear h; 
    clear time;
    clear yFOR;

end %Ends looping over different h-values


%
%Plots CONVERGENCE RESULTS in PAPER
%
plots_For_Convergence_Paper(Y_SOL,Y_SOL2,NVecPlot,hVecPlot,tE,tS,y0,errMat,errMat2,hVec,err_h,err_h2,NVec,timeV);


fprintf('\nWelp, thats it folks!\n\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Makes plots for paper
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plots_For_Convergence_Paper(Y_Sol,Y_Sol2,NVec2,hVec2,tE,tS,y0,errMat,errMat2,hVec,err_h,err_h2,Nvec,time)

%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Attributes
%%%%%%%%%%%%%%%%%%%%%%%%%
fs = 18; % FontSize
ms = 35; % MarkerSize
lw = 6;  % LineWidth

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots the EXACT vs. NUMERICAL Solution for various time-steps, h
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)

%
% Colors for plotting
strColor = [0 0 0; 
    0 0.2 0.7;
    0 0.447 0.9;
    0.3010 0.7450 0.9330;
    0.5940, 0.0940, 0.6560];

for i=1:length(Y_Sol2(1,:))
    h = hVec2(i);
    hVecStr{i} = strcat('\Deltat = ',num2str(h));
    N=  NVec2(i);
    t = tS:h:tE;
    %str = strColor{i};
    str = strColor(i,:);
    plot(t,Y_Sol2(1:N+1,i),'-','LineWidth',lw,'Color',str); hold on;
end
t=tS:0.005:tE;
for i=1:length(t)
    yExact(i) = Exact(t(i),y0,2);
end
plot(t,yExact,'r-','LineWidth',lw); hold on;
leg=legend(hVecStr{1},hVecStr{2},hVecStr{3},hVecStr{4},hVecStr{5},'Exact Solution','Location','NorthEast');
%title('Exact vs Numerical Solutions');
xlabel('t');
ylabel('y(t)');
set(gca,'FontSize',fs);
set(leg,'FontSize',fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots Convergence Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
loglog(hVec,err_h2,'b.','MarkerSize',ms); hold on;
xlabel('\Deltat');
ylabel('Absolute Error');
%title('Convergence Study: y(t) = sin(2*pi*t)+1');
set(gca,'FontSize',fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots Computational Time Study
%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
loglog(Nvec,time,'r.-','LineWidth',lw,'MarkerSize',ms); hold on;
xlabel('# of Time-Steps btwn [tS,tE]/\Deltat');
ylabel('Computational Time');
%title('Computational Time Study');
set(gca,'FontSize',fs);
set(leg,'FontSize',fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots Computational Time Study
%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
loglog(2./Nvec,time,'r.-','LineWidth',lw,'MarkerSize',ms); hold on;
xlabel('\Deltat');
ylabel('Computational Time');
%title('Computational Time Study');
set(gca,'FontSize',fs);
set(leg,'FontSize',fs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function that computes Error
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function error = compute_Error(time,y,y0,flag)

%Computes Error During Simulation

exact_sol = zeros(1,length(time));
for i=1:length(time)
   exact_sol(i) = Exact(time(i),y0,flag); 
end

error = y - exact_sol;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the RHS of the ODE: y' = f(t,y)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = f(t,y,flag)

%y' = f(t,y)

if flag == 1
    val = y;
elseif flag==2
    val = 2*pi*cos(2*pi*t);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the Exact Sol'n
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = Exact(t,y0,flag)

%Exact sol'n to ODE
if flag == 1
    val = y0*exp(t);
elseif flag == 2
    val = sin(2*pi*t)+1; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function tells you how many time-steps you have for the simulation
% (This is used for the convergence study)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NVec = give_Me_NumberOfGridPts()

N1 = 1:1:9;
N2 = 10:10:90;
N3 = 100:100:900;
N4 = 1000:1000:9*1e3;
N5 = 1e4:1e4:9*1e4;
N6 = 1e5:1e5:9*1e5;
N7 = 1e6:1e6:9*1e6;
NVec = [N1 N2 N3 N4 N5 N6 N7];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function that prints the info about the simulation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_info()

fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');

fprintf('\n Author: Nicholas A. Battista\n');
fprintf(' Created: August 12, 2014\n');
fprintf(' Date of Last Revision: February 3, 2019\n\n');

fprintf(' This function solves the following ODE:\n\n');
fprintf(' dy/dt = f(t,y)\n');
fprintf(' y(0) = y0\n\n');
fprintf(' using Eulers Method. It then does a convergence study for various h values \n\n');
fprintf(' Note: It performs the convergence study for the known ODEs,\n\n');

fprintf(' dy/dt = y, with y(0) = 1; so exact solution is: y(t) = e^t\n\n');
fprintf(' and\n\n');
fprintf(' dy/dt = 2*pi*cos(2*pi*t), w/ y(0)=1 w/ exact solution y(t) = sin(2*pi*t)+1\n\n');

fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');

