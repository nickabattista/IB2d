%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  battistn@tcnj.edu
% Date Created: May 27th, 2015
% Date Modified: February 15, 2020
% Institution Created: UNC-CH
% Institution Modified: TCNJ
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
% If you would like us to add a specific muscle model, please let Nick (battistn[@]tcnj.edu) know.
%
%--------------------------------------------------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: updates the target point positions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function targets = update_Target_Point_Positions(dt,current_time,targets)

% last index of moving target points
last_id=284;

IDs = targets(1:last_id,1);         % Stores Lag-Pt IDs in col vector
xPts= targets(1:last_id,2);         % Original x-Values of x-Target Pts.
yPts= targets(1:last_id,3);         % Original y-Values of y-Target Pts.
%kStiffs = targets(:,4);            % Stores Target Stiffnesses 



% Coefficients for Polynomial Phase-Interpolation
a = 2.739726027397260;  % y1(t) = at^2
b = 2.739726027397260;  % y3(t) = -b(t-1)^2+1
c = -2.029426686960933; % y2(t) = ct^3 + dt^2 + gt + h
d = 3.044140030441400;
g = -0.015220700152207;
h = 0.000253678335870;

% Period Info
tP1 = 0.05;  % Interpolation: P1->P2 
tP2 = 0.01;  % Waiting in P2 Phase
tP3 = 0.05;  % Interpolation: P2->P3
tP4 = 0.01;  % Waiting in P3 Phase

period = tP1+tP2+tP3+tP4;          % Period
t = rem(current_time,period);      % Current time in simulation ( 'modular arithmetic to get time in period')


% Read In Pts!
[xP1,yP1,xP2,yP2,xP3,yP3,~,~] = read_File_In('All_Positions.txt',last_id);


 
if (t <= tP1) 

        %PHASE 1 --> PHASE 2

        t1 = 0.1*tP1;   
        t2 = 0.9*tP1;
        if (t<=t1) 							%For Polynomial Phase Interp.
            g1 = a*power((t/tP1),2);
        elseif ((t>=t1)&&(t<t2)) 
            g1 = c*power((t/tP1),3) + d*power((t/tP1),2) + g*(t/tP1) + h;
        elseif (t>=t2)
            g1 = -b*power(((t/tP1) - 1),2) + 1;
        end

        xPts(IDs) = xP1 + g1*( xP2 - xP1 );
        yPts(IDs) = yP1 + g1*( yP2 - yP1 );	

        %if ( ( t > tP1-dt/2) && ( t < tP1+dt/2) )
        %    xM1 = xPts(IDs);
        %    yM1 = yPts(IDs);
        %    please_Print_Vertices_To_File(xP1,yP1,xP2,yP2,xP3,yP3,xM1,yM1);
        %end

elseif ((t>tP1)&&(t<=(tP1+tP2)))

        % PHASE 2's * Waiting Period *

        tprev = tP1;
        t1 = 0.1*tP2 + tprev;
        t2 = 0.9*tP2 + tprev;
        if (t<t1) 							%//For Polynomial Phase Interp.
            g2 = a*power( ( (t-tprev)/tP2) ,2);
        elseif ((t>=t1)&&(t<t2)) 
            g2 = c*power( ( (t-tprev)/tP2) ,3) + d*power( ((t-tprev)/tP2) ,2) + g*( (t-tprev)/tP2) + h;
        elseif (t>=t2) 
            g2 = -b*power( (( (t-tprev)/tP2) - 1) ,2) + 1;
        end			

        xPts(IDs) = xP2 + g2*( xP2 - xP2 );
        yPts(IDs) = yP2 + g2*( yP2 - yP2 );            

        %if (( t > tP1+tP2-dt) && ( t < tP1+tP2+dt) )
        %    xM1 = xPts(IDs);
        %    yM1 = yPts(IDs);
        %    please_Print_Vertices_To_File(xP1,yP1,xP2,yP2,xP3,yP3,xM1,yM1);
        %end



elseif ((t>tP1+tP2)&&(t<=(tP1+tP2+tP3)))

        %PHASE 2 --> PHASE 3 

        tprev = tP1+tP2;
        t1 = 0.1*tP3 + tprev;
        t2 = 0.9*tP3 + tprev;
        if (t<t1) 							%//For Polynomial Phase Interp.
            g2 = a*power( ( (t-tprev)/tP3) ,2);
        elseif ((t>=t1)&&(t<t2)) 
            g2 = c*power( ( (t-tprev)/tP3) ,3) + d*power( ((t-tprev)/tP3) ,2) + g*( (t-tprev)/tP3) + h;
        elseif (t>=t2) 
            g2 = -b*power( (( (t-tprev)/tP3) - 1) ,2) + 1;
        end			

        xPts(IDs) = xP2 + g2*( xP3 - xP2 );
        yPts(IDs) = yP2 + g2*( yP3 - yP2 );

        %if (( t > tP1+tP2+tP3-dt ) && ( t < tP1+tP2+tP3+dt ) )
        %    xM1 = xPts(IDs);
        %    yM1 = yPts(IDs);
        %    please_Print_Vertices_To_File(xP1,yP1,xP2,yP2,xP3,yP3,xM1,yM1);
        %end

elseif ((t>tP1+tP2+tP3)&&(t<=(tP1+tP2+tP3+tP4)))

        % PHASE 3's * Waiting Period *

        tprev = tP1+tP2+tP3;
        t1 = 0.1*tP4 + tprev;
        t2 = 0.9*tP4 + tprev;
        if (t<t1) 							%//For Polynomial Phase Interp.
            g2 = a*power( ( (t-tprev)/tP4) ,2);
        elseif ((t>=t1)&&(t<t2)) 
            g2 = c*power( ( (t-tprev)/tP4) ,3) + d*power( ((t-tprev)/tP4) ,2) + g*( (t-tprev)/tP4) + h;
        elseif (t>=t2) 
            g2 = -b*power( (( (t-tprev)/tP4) - 1) ,2) + 1;
        end			

        xPts(IDs) = xP3 + g2*( xP3 - xP3 );
        yPts(IDs) = yP3 + g2*( yP3 - yP3 );

        %if (( t > tP1+tP2+tP3+tP4-dt) && ( t < tP1+tP2+tP3+tP4+dt) )
        %    xM1 = xPts(IDs);
        %    yM1 = yPts(IDs);
        %    please_Print_Vertices_To_File(xP1,yP1,xP2,yP2,xP3,yP3,xM1,yM1);
        %end

end

targets(IDs,2) = xPts(IDs); % Store new xVals
targets(IDs,3) = yPts(IDs); % Store new yVals



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints all Vertices to File
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function please_Print_Vertices_To_File(X1,Y1,X2,Y2,X3,Y3,xM,yM)

delete('All_Positions.txt');

fileID = fopen('All_Positions.txt','w');
for j=1:length(X1)
    fprintf(fileID,'%1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e\n', X1(j),Y1(j),X2(j),Y2(j),X3(j),Y3(j),xM(j),yM(j) );
end
fclose(fileID);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: reads in info from file that contains probability distributions
% (rows) for each game of Bingo w/ N players (columns)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x1,y1,x2,y2,x3,y3,xM,yM] = read_File_In(file_name,last_id)

filename = file_name;  %Name of file to read in

fileID = fopen(filename);

    % Read in the file, use 'CollectOutput' to gather all similar data together
    % and 'CommentStyle' to to end and be able to skip lines in file.
    C = textscan(fileID,'%f %f %f %f %f %f %f %f','CollectOutput',1);

fclose(fileID);        %Close the data file.

mat_info = C{1};   %Stores all read in data

%Store all elements in matrix
mat = mat_info(1:end,1:end);

x1 =  mat(1:last_id,1); % store x-Phase 1 positions
y1 =  mat(1:last_id,2); % store y-Phase 1 positions
x2 =  mat(1:last_id,3); % store x-Phase 2 positions
y2 =  mat(1:last_id,4); % store y-Phase 2 positions
x3 =  mat(1:last_id,5); % store x-Phase 3 positions
y3 =  mat(1:last_id,6); % store y-Phase 3 positions
xM =  mat(1:last_id,7); % store x-Phase 3 positions
yM =  mat(1:last_id,8); % store y-Phase 3 positions