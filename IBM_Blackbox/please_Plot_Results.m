%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Plots the Lagrangian structure with:
%           (1): the background velocity field
%           (2): the Lagrangian structure itself
%           (3): the vorticity field in a colormap
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function please_Plot_Results(dx,dy,X,Y,U,V,chiX,chiY)

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

diffX = abs([chiX(2:end);chiX(1)]-chiX);
diffY = abs([chiY(2:end);chiY(1)]-chiY);

locX = find(diffX > Lx/2); % assuming x value can't change more than Lx/2 (without crossing boundary)
locY = find(diffY > Ly/2); % assuming y value can't change more than Ly/2 (without crossing boundary)
loc =  sort(unique([locX;locY]));

clf;


% % % % PLOTS THE VELOCITY (EULERIAN) GRID AND LAGRANGIAN MESH % % % % 

subplot(1,3,1)
axis([0 Lx 0 Ly]);
title('VELOCITY');
xlabel('x'); ylabel('y');
hold all;

quiver(X,Y,U,V);

   loc = [0;loc;length(chiX)];
   for i=2:length(loc)
       plot(chiX(loc(i-1)+1:loc(i)),chiY(loc(i-1)+1:loc(i)),'m','LineWidth',3);
   end

axis square;


 
% % % % % % % % PLOTS THE LAGRANGIAN MESH! % % % % % % % %


subplot(1,3,2)
axis([0 Lx 0 Ly]);
title('LAGRANGIAN PTS');
xlabel('x'); ylabel('y');
hold all;

   loc = [0;loc;length(chiX)];
   for i=2:length(loc)
       plot(chiX(loc(i-1)+1:loc(i)),chiY(loc(i-1)+1:loc(i)),'m','LineWidth',3);
   end

axis square;



% % % % % PLOTS THE VORTICITY AND LAGRANGIAN MESH! % % % % %

subplot(1,3,3)
axis([0 Lx 0 Ly]);
title('VORTICITY');
xlabel('x'); ylabel('y'); 
hold all;

%Compute Vorticity and Plot It against Lagrangian Grid!
vort = give_Me_Vorticity(U,V,dx,dy); 
vort = vort';
x = X(1,:); y = Y(:,1);
contourf(x,y,flipud(rot90(vort)),10); hold on;

   loc = [0;loc;length(chiX)];
   for i=2:length(loc)
       plot(chiX(loc(i-1)+1:loc(i)),chiY(loc(i-1)+1:loc(i)),'m','LineWidth',3);
   end


axis square;
drawnow;

hold off;
set(gca,'Box','on');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Computes vorticity from two matrices, U and V, where each
% matrix is the discretized field of velocity values either for x or y,
% respectively.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vort = give_Me_Vorticity(U,V,dx,dy)

% w = ( dv/dx - du/dy )\hat{z}

%Compute dv/dx using central differencing! (maintains periodicity)
dvdx = D(V,dx,'x');

%Compute du/dy using central differencing! (maintains periodicity)
dudy = D(U,dy,'y');
    
vort = ( dvdx - dudy );