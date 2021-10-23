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
%	4. Muscle-Model (combined Force-Length-Velocity model, "Hill+(Length-Tension)")
%
% One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting lengths, etc
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
% 
% If you would like us to add a specific muscle model, please let Nick (nick.battista@unc.edu) know.
%
%--------------------------------------------------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: updates the target point positions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function targets = update_Target_Point_Positions(dt,current_time,targets)


IDs = targets(:,1);                 % Stores Lag-Pt IDs in col vector
%xPts= targets(:,2);                % Original x-Values of x-Target Pts.
yPts= targets(:,3);                 % Original y-Values of y-Target Pts.
%kStiffs = targets(:,4);             % Stores Target Stiffnesses 

N_target = length(targets(:,1));    %Gives total number of target pts!

% NEEDED FOR UPDATE FILE FROM CORAL GEOMETRY! %
N_half = 45;                        % # of pts. on ONE tentacle

    %frames per second in data
     F = 30.0; 
     %pixels per dm in data 
     S = 1733.333333333333333333;
     %time in frames (mod one pulse)
     t =  F*mod(current_time,2.30);
     tp = F*2.25; 
     tend = F*2.30; 
     %initial offset from center of the domain in dm
     xoffset = 0.019; 
     yoffset = 0.25; 
     xoffset2grid = 0.3;

     qc11 = -8.018203211249667e-13;
     qc12 = 1.706310971452060e-10;
     qc13 = -1.371526883228119e-8;
     qc14 = 5.116506115621489e-7;
     qc15 = -8.613944795521634e-6;
     qc16 = 4.951095471423292e-5;
     qc17 = -1.237254488643602e-4;
     qc21 = -1.093988700839764e-10;
     qc22 = 2.370737151797367e-8;
     qc23 = -1.946684756263760e-6;
     qc24 = 7.447999061932255e-5;
     qc25 = -1.285757196003425e-3;
     qc26 = 7.184871148722453e-3;
     qc27 = -2.357284935085283e-2;    
     qc31 = -3.111491615328405e-9;
     qc32 = 6.922752279229525e-7; 
     qc33 = -5.951659336196739e-5; 
     qc34 = 2.461541908755476e-3; 
     qc35 = -4.830999119980220e-2; 
     qc36 = 3.293063238297988e-1;
     qc37 = -1.893567487439901;
     qc41 = 2.028150369106308e-8;
     qc42 = -4.188949936239104e-6; 
     qc43 = 3.237380275082298e-4;
     qc44 = -1.158763275070524e-2; 
     qc45 = 1.908244836752949e-1;
     qc46 = -1.142372536532416; 
     qc47 = 3.584246435337805;

     
    if (t<=tp)

        qpoly1=qc17+qc16*t+qc15*t*t+qc14*t*t*t+qc13*t*t*t*t+qc12*t*t*t*t*t+qc11*t*t*t*t*t*t; 
        qpoly2=qc27+qc26*t+qc25*t*t+qc24*t*t*t+qc23*t*t*t*t+qc22*t*t*t*t*t+qc21*t*t*t*t*t*t;
        qpoly3=qc37+qc36*t+qc35*t*t+qc34*t*t*t+qc33*t*t*t*t+qc32*t*t*t*t*t+qc31*t*t*t*t*t*t;
        qpoly4=qc47+qc46*t+qc45*t*t+qc44*t*t*t+qc43*t*t*t*t+qc42*t*t*t*t*t+qc41*t*t*t*t*t*t;
	
    else 
        qpolytf1=qc17+qc16*tp+qc15*tp*tp+qc14*tp*tp*tp+qc13*tp*tp*tp*tp+qc12*tp*tp*tp*tp*tp+qc11*tp*tp*tp*tp*tp*tp;
        qpolytf2=qc27+qc26*tp+qc25*tp*tp+qc24*tp*tp*tp+qc23*tp*tp*tp*tp+qc22*tp*tp*tp*tp*tp+qc21*tp*tp*tp*tp*tp*tp;
        qpolytf3=qc37+qc36*tp+qc35*tp*tp+qc34*tp*tp*tp+qc33*tp*tp*tp*tp+qc32*tp*tp*tp*tp*tp+qc31*tp*tp*tp*tp*tp*tp;
        qpolytf4=qc47+qc46*tp+qc45*tp*tp+qc44*tp*tp*tp+qc43*tp*tp*tp*tp+qc42*tp*tp*tp*tp*tp+qc41*tp*tp*tp*tp*tp*tp;

        qpoly1 = qpolytf1*((tend-tp)-(t-tp))/(tend-tp)+qc17*(t-tp)/(tend-tp); 
        qpoly2 = qpolytf2*((tend-tp)-(t-tp))/(tend-tp)+qc27*(t-tp)/(tend-tp);
        qpoly3 = qpolytf3*((tend-tp)-(t-tp))/(tend-tp)+qc37*(t-tp)/(tend-tp);
        qpoly4 = qpolytf4*((tend-tp)-(t-tp))/(tend-tp)+qc47*(t-tp)/(tend-tp);
    end
     
     
     
for i=1:N_target                    % Loops over all target points!
    
    shifted_Xt1 = yPts(i)-yoffset;
   
    X_target = xoffset+qpoly4/S-qpoly3*shifted_Xt1+qpoly2*S*shifted_Xt1*shifted_Xt1-qpoly1*S*S*shifted_Xt1*shifted_Xt1*shifted_Xt1;

    if i<=N_half
        X_target = -X_target+ xoffset2grid;
    else
        X_target = X_target + xoffset2grid;
    end
    
    targets(IDs(i),2) = X_target;      % Store new xVals
    %targets(IDs(i),3) = yPts(IDs(i)); % Store new yVals

end

