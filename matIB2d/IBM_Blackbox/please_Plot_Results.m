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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Plots the Lagrangian structure with:
%           (1): the background velocity field
%           (2): the Lagrangian structure itself
%           (3): the vorticity field in a colormap
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [loc, diffy] =  please_Plot_Results(ds,X,Y,U,V,vort,uMag,p,C,chiX,chiY,lagPlot,velPlot,vortPlot,pressPlot,uMagPlot,conPlot,firstPrint,loc,diffy,spacing)

%X,Y:  (x,y) values
%U,V:  x-directed, y-directed velocities respectively
%vort: vorticity
%uMag: magnitude of velocity
%p:    pressure

%FLAGS FOR PLOTTING:
%pressPlot - if you want pressure plot = 1
%uMagPlot  - if you want mag. velocity plot = 1
%vortPlot  - if you want vorticity plot = 1
%velPlot   - if you want velocity plot = 1
%lagPlot   - if you want lag. point ONLY plot = 1
%conPlot   - if you want concentration plot =1
%
% Assumption: Assuming chiX and chiY are column vectors
% Assumption: Assuming chiX(i+1)-chiX(i) < .5 and chiY(i+1)-chiY(i) < .5, for all points that don't cross the boundary
%
% This script has been adapted from the Mat-IB code by Hiens and Stockie, 2014.
%

Lx = X(1,end)+X(1,2);
Ly = Y(end,1)+Y(2,1);

%Shift so inside of interval [0,Lx] or [0,Ly]
chiX = mod(chiX,Lx);
chiY = mod(chiY,Ly);

%Connect Geometry Appropriately At Beginning Of Simulation then Store it!
if firstPrint == 1
    diffX = abs([chiX(2:end);chiX(1)]-chiX);
    diffY = abs([chiY(2:end);chiY(1)]-chiY);

    locX = find(diffX > spacing ); % assuming x value can't change more than Lx/2 (without crossing boundary)
    locY = find(diffY > spacing ); % assuming y value can't change more than Ly/2 (without crossing boundary)
    loc =  sort(unique([locX;locY]));

    diffy = sqrt( (chiX(end)-chiX(1) )^2 + ( chiY(end)-chiY(1) )^2 );
end

clf; %Clear previous plots :)




figure(1) 
numPlots = lagPlot+velPlot+vortPlot+pressPlot+uMagPlot+conPlot;

ct = 1;

% % % % % PLOTS LAGRANGIAN POINTS ONLY (if selected) % % % % %

if lagPlot == 1
    subplot(1,numPlots,ct)
    axis([0 Lx 0 Ly]);
    title('LAGRANGIAN PTS');
    xlabel('x'); ylabel('y');
    hold all;

    loc = [0;loc;length(chiX)];
    for i=2:length(loc)
       plot(chiX(loc(i-1)+1:loc(i)),chiY(loc(i-1)+1:loc(i)),'m','LineWidth',3);
    end
    if diffy < 5*ds
        xTemp = [chiX(1) chiX(end)];
        yTemp = [chiY(1) chiY(end)];
        plot(xTemp(1:2),yTemp(1:2),'m','LineWidth',3);
    end

    axis square;

    ct=ct+1;
end

% % % % % PLOTS VORTICITY + LAG Pts. (if selected) % % % % %

if vortPlot == 1
    
    subplot(1,numPlots,ct)
    %
    axis([0 Lx 0 Ly]);
    title('VORTICITY');
    xlabel('x'); ylabel('y'); 
    hold all;

    %Compute Vorticity and Plot It against Lagrangian Grid!
    x = X(1,:); y = Y(:,1);
    contourf(x,y,flipud(rot90(vort)),10); hold on;

    loc = [0;loc;length(chiX)];
    for i=2:length(loc)
        plot(chiX(loc(i-1)+1:loc(i)),chiY(loc(i-1)+1:loc(i)),'m','LineWidth',3);
    end
    if diffy < 5*ds
        xTemp = [chiX(1) chiX(end)];
        yTemp = [chiY(1) chiY(end)];
        plot(xTemp(1:2),yTemp(1:2),'m','LineWidth',3);
    end
    
    axis square;

    ct=ct+1;
end

% % % % % PLOTS PRESSURE + LAG Pts. (if selected) % % % % %

if pressPlot == 1
    
    subplot(1,numPlots,ct)
    %
    axis([0 Lx 0 Ly]);
    title('PRESSURE');
    xlabel('x'); ylabel('y'); 
    hold all;

    %Use Pressure and Plot It against Lagrangian Grid!
    x = X(1,:); y = Y(:,1);
    contourf(x,y,p,6); hold on;

    loc = [0;loc;length(chiX)];
    for i=2:length(loc)
        plot(chiX(loc(i-1)+1:loc(i)),chiY(loc(i-1)+1:loc(i)),'m','LineWidth',3);
    end
    if diffy < 5*ds
        xTemp = [chiX(1) chiX(end)];
        yTemp = [chiY(1) chiY(end)];
        plot(xTemp(1:2),yTemp(1:2),'m','LineWidth',3);
    end

    axis square;
    
    ct=ct+1;
end

% % % % % PLOTS MAGNITUDE OF VELOCITY + LAG Pts. (if selected) % % % % %

if uMagPlot == 1
    
    subplot(1,numPlots,ct)
    %
    axis([0 Lx 0 Ly]);
    title('MAGNITUDE OF VELOCITY');
    xlabel('x'); ylabel('y'); 
    hold all;

    %Use Mag. Velocity and Plot It against Lagrangian Grid!
    x = X(1,:); y = Y(:,1);
    contourf(x,y,uMag,6); hold on;

    loc = [0;loc;length(chiX)];
    for i=2:length(loc)
        plot(chiX(loc(i-1)+1:loc(i)),chiY(loc(i-1)+1:loc(i)),'m','LineWidth',3);
    end
    if diffy < 5*ds
        xTemp = [chiX(1) chiX(end)];
        yTemp = [chiY(1) chiY(end)];
        plot(xTemp(1:2),yTemp(1:2),'m','LineWidth',3);
    end

    axis square;
    
    ct=ct+1;
end

% % % % % PLOTS VELOCITY FIELD + LAG Pts. (if selected) % % % % %

if velPlot == 1
    
    subplot(1,numPlots,ct)
    %
    axis([0 Lx 0 Ly]);
    title('VELOCITY');
    xlabel('x'); ylabel('y');
    hold all;
    
    
    quiver(X,Y,U,V); %Print Velocity Field
    
    
    loc = [0;loc;length(chiX)];
    for i=2:length(loc)
       plot(chiX(loc(i-1)+1:loc(i)),chiY(loc(i-1)+1:loc(i)),'m','LineWidth',3);
    end
    if diffy < 5*ds
        xTemp = [chiX(1) chiX(end)];
        yTemp = [chiY(1) chiY(end)];
        plot(xTemp(1:2),yTemp(1:2),'m','LineWidth',3);
    end

    axis square;
    
    ct=ct+1;
end

if conPlot == 1
    
    subplot(1,numPlots,ct)
    %
    axis([0 Lx 0 Ly]);
    title('Concentration');
    xlabel('x'); ylabel('y'); 
    hold all;

    %Compute Vorticity and Plot It against Lagrangian Grid!
    x = X(1,:); y = Y(:,1);
    %contourf(x,y,flipud(rot90(vort)),10); hold on;
    pcolor(X,Y,C);
    shading interp;

    loc = [0;loc;length(chiX)];
    for i=2:length(loc)
        plot(chiX(loc(i-1)+1:loc(i)),chiY(loc(i-1)+1:loc(i)),'m','LineWidth',3);
    end
    if diffy < 5*ds
        xTemp = [chiX(1) chiX(end)];
        yTemp = [chiY(1) chiY(end)];
        plot(xTemp(1:2),yTemp(1:2),'m','LineWidth',3);
    end
    
    axis square;

    ct=ct+1;
end

drawnow;

hold off;
set(gca,'Box','on');
