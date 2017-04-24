%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Runs Advection-Diffusion Code
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function adv_diff()

L = 1;              % Size of computational grid
N = 64;            % # of grid points
dx = L/(N-1);dy=dx; % Grid Resolution
k = 1e-6;           % Diffusive Coefficient

dt = 1e-3;          % time-step
Tfinal = 40.0;      % final time
tVec = 0:dt:Tfinal; % time-vector

[uX,uY] = give_Velocity_Fields(N,L,dx); % Give background velocity fields
plot_Vector_Field(L,uX,uY,dx);
C = give_Initial_Concentration(N,L,dx); % Give initial concentration

% SAVING DATA TO VTK %
print_dump = 100;
ctsave = 0;
% CREATE VIZ_IB2D FOLDER and VISIT FILES
mkdir('viz_IB2d');
print_vtk_files(ctsave,uX',uY',L,N,C)

fprintf('\nStarting the advection-diffusion time-stepping!\n');

%
% START THE TIME-STEPPING
%
for i=1:length(tVec)

    % Update the advection-diffusion
    
    % Flux Limiter
    C = please_Update_Adv_Diff_Concentration_Flux_Limiter(C,dt,dx,dy,uX,uY,k);
    
    % Split (temporally) + Upwind
    %C = please_Update_Adv_Diff_Concentration_Split(C,dt,dx,dy,uX,uY,k);
    
    % Split (temporally) + Lax-Wendroff
    %C = please_Update_Adv_Diff_Concentration_Split_Lax_Wendroff(C,dt,dx,dy,uX,uY,k);

    % Un-Split (temporally) + Upwind
    %C = please_Update_Adv_Diff_Concentration_Unsplit(C,dt,dx,dy,uX,uY,k);


    
    % Save files info!
    ctsave = ctsave + 1;
    if mod(ctsave,print_dump) == 0
        print_vtk_files(ctsave,uX',uY',L,N,C);
    end
    fprintf('Time: %d\n',tVec(i));
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: give initial concentration!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C = give_Initial_Concentration(N,L,dx)

%WHERE OUTER TUBE LIES
%xMin = 0.15; xMax = 0.45;
%yMin = 0.85; yMax = 1.15;

xMin = 0.6; xMax = 0.8;
yMin = 0.2; yMax = 0.4;

xMid = (xMin+xMax)/2;
yMid = (yMin+yMax)/2;

xDiff = (xMax-xMin)/2;
yDiff = (yMax-yMin)/2;

x = 0:dx:L;
y = 0:dx:L;

inds = give_Me_Indices_To_Apply_Force(x,y,xMin,xMax,yMin,yMax);

C = zeros(N,N);
for i=1:length( inds(:,1) )
    xInd = inds(i,1);
    yInd = inds(i,2);
    xPt = x(xInd);
    yPt = y(yInd);
    C(yInd,xInd ) = (-0.5/yDiff^2)*( (yPt-yMid) - yDiff )*( (yPt-yMid) + yDiff ) +  (-0.5/xDiff^2)*( (xPt-xMid) - xDiff )*( (xPt-xMid) + xDiff ); %1.0;
    %C(xInd,yInd ) = (-1.0/xDiff^2)*( (xPt-xMid) - xDiff )*( (xPt-xMid) + xDiff ); %1.0;

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: give prescribed velocity fields!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [uX,uY] = give_Velocity_Fields(N,L,dx)

xGrid = 0:dx:L; yGrid=xGrid; % Construct grid values
xC = L/2; yC = L/2;          % Center point

uX = zeros(N,N); uY = uX;    % Allocate memory

for i=1:length(xGrid)
    
    x = xGrid(i);
    
    for j=1:length(yGrid)
    
        y = yGrid(j);
        
        uX(j,i) = -( y - yC )/5;
        uY(j,i) = ( x - xC )/5;
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: plot vector field
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_Vector_Field(L,uX,uY,dx)

figure(55)
xGrid = 0:dx:L; yGrid=xGrid; % Construct grid values
%
axis([0 L 0 L]);
title('VELOCITY');
xlabel('x'); ylabel('y');
hold all;

quiver(xGrid,yGrid,uX,uY); %Print Velocity Field

drawnow;
    