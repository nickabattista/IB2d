%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
% 	fluid-structure eraction models. This version of the code is based off of
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
% 	3. Target Pos
%	4. Muscle-Model (combined Force-Length-Velocity model, "Hill+(Length-Tension)")
%
% One is able to update those Lagrangian Structure parameters, e.g., spring ants, resting lengths, etc
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
% 
% If you would like us to add a specific muscle model, please let Nick (nick.battista@unc.edu) know.
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: updates the beam attributes!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function beams_info = update_nonInv_Beams(dt,current_time,beams_info)


% beams_info:   col 1: 1ST PT.
%               col 2: MIDDLE PT. (where force is exerted)
%               col 3: 3RD PT.
%               col 4: beam stiffness
%               col 5: x-curvature
%               col 6: y-curvature

%IDs = beams_info(:,1);   % Gives Middle Pt.

% Coefficients for Polynomial Phase-Interpolation
a = 2.739726027397260;  % y1(t) = at^2
b = 2.739726027397260;  % y3(t) = -b(t-1)^2+1
c = -2.029426686960933; % y2(t) = ct^3 + dt^2 + gt + h
d = 3.044140030441400;
g = -0.015220700152207;
h = 0.000253678335870;

% Period Info
tP1 = 0.25;                       % Down-stroke
tP2 = 0.25;                       % Up-stroke
period = tP1+tP2;                  % Period
t = rem(current_time,period);      % Current time in simulation ( 'modular arithmetic to get time in period')

% Read In y_Pts for two Phases!
[xP1,yP1,yP2] = read_File_In('swimmer.phases'); % NOTE xP1 = xP2 
xP2 = xP1;

%
% FIRST WE COMPUTE THE INTERPOLATE GEOMETRY BETWEEN BOTH PHASES
%

    %PHASE 1 --> PHASE 2
    if (t <= tP1) 

        %tprev = 0.0;
        t1 = 0.1*tP1;   
        t2 = 0.9*tP1;
        if (t<t1) 							%For Polynomial Phase Interp.
            g1 = a*power((t/tP1),2);
        elseif ((t>=t1)&&(t<t2)) 
            g1 = c*power((t/tP1),3) + d*power((t/tP1),2) + g*(t/tP1) + h;
        elseif (t>=t2)
            g1 = -b*power(((t/tP1) - 1),2) + 1;
        end
			
        xPts = xP1 + g1*( xP2 - xP1 );	
        yPts = yP1 + g1*( yP2 - yP1 );	
		
    %PHASE 2 --> PHASE 1
    elseif ((t>tP1)&&(t<=(tP1+tP2)))
			
        tprev = tP1;
        t1 = 0.1*tP2 + tP1;
        t2 = 0.9*tP2 + tP1;
        if (t<t1) 							%//For Polynomial Phase Interp.
            g2 = a*power( ( (t-tprev)/tP2) ,2);
        elseif ((t>=t1)&&(t<t2)) 
            g2 = c*power( ( (t-tprev)/tP2) ,3) + d*power( ((t-tprev)/tP2) ,2) + g*( (t-tprev)/tP2) + h;
        elseif (t>=t2) 
            g2 = -b*power( (( (t-tprev)/tP2) - 1) ,2) + 1;
        end			

        xPts = xP2 + g2*( xP1 - xP2 );
        yPts = yP2 + g2*( yP1 - yP2 );
    
    end


%
% NOW WE UPDATE THE CURAVTURES APPROPRIATELY
%
beams_info(:,5) = xPts(1:end-2)+xPts(3:end)-2*xPts(2:end-1);
beams_info(:,6) = yPts(1:end-2)+yPts(3:end)-2*yPts(2:end-1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: reads in info from file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x1,y1,y2] = read_File_In(file_name)

filename = file_name;  %Name of file to read in

fileID = fopen(filename);

    % Read in the file, use 'CollectOutput' to gather all similar data together
    % and 'CommentStyle' to to end and be able to skip lines in file.
    C = textscan(fileID,'%f %f %f','CollectOutput',1);

fclose(fileID);        %Close the data file.

mat_info = C{1};   %Stores all read in data

%Store all elements in matrix
mat = mat_info(1:end,1:end);

x1 =  mat(:,1); %store xVals1/2
y1 =  mat(:,2); %store yVals1 
y2 =  mat(:,3); %store yVals2