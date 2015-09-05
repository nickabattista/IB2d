%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Setting up advection-diffusion solver
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C = please_Update_Adv_Diff_Concentration(C,dt,dx,dy,uX,uY,k)

% C:     concentration 
% dt:    time-step
% dx,dy: spatial steps in x and y, respectively
% uX:    x-Component of Velocity
% uY:    y-Component of Velocity
% k:     diffusion coefficient

% Compute Necessary Derivatives 
Cx = D(C,dx,'x');
Cy = D(C,dy,'y');
Cxx = DD(C,dx,'x');
Cyy = DD(C,dy,'y');
    
% Update Concentration 
% C = C + dt * ( k*(Cxx+Cyy) - uX'.*Cx - uY'.*Cy );

C = C + dt * ( k*(Cxx+Cyy) - uX.*Cx - uY.*Cy );

