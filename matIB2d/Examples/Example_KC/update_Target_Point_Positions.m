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
xPts= targets(:,2);                 % Original x-Values of x-Target Pts.
yPts= targets(:,3);                 % Original y-Values of y-Target Pts.
kStiffs = targets(:,4);             % Stores Target Stiffnesses 

N_target = length(targets(:,1));    %Gives total number of target pts!


% Coefficients for Polynomial Phase-Interpolation
a = 2.739726027397260;  % y1(t) = at^2
b = 2.739726027397260;  % y3(t) = -b(t-1)^2+1
c = -2.029426686960933; % y2(t) = ct^3 + dt^2 + gt + h
d = 3.044140030441400;
g = -0.015220700152207;
h = 0.000253678335870;

% Period Info
tP1 = 0.05;                        % Contraction time
tP2 = 0.05;                       % Expansion time
tP3 = 0.05;
tP4 = 0.05;
tP5 = 0.05;
period = tP1+tP2+tP3+tP4;           % Period
t = rem(current_time,period);       % Current time in simulation ( 'modular arithmetic to get time in period')


% Read In Pts!
[xP1,yP1,xP2,yP2,xP3,yP3,xM,yM] = read_File_In('All_Positions.txt');



for i=1:N_target                    % Loops over all target points!
    
    
    if (t <= tP1) 

			%PHASE 1 --> PHASE 2
			
			tprev = 0.0;
			t1 = 0.1*tP1;   
			t2 = 0.9*tP1;
			if (t<=t1) 							%For Polynomial Phase Interp.
				g1 = a*power((t/tP1),2);
            elseif ((t>=t1)&&(t<t2)) 
				g1 = c*power((t/tP1),3) + d*power((t/tP1),2) + g*(t/tP1) + h;
            elseif (t>=t2)
				g1 = -b*power(((t/tP1) - 1),2) + 1;
            end
			
			xPts(IDs(i)) = xP1(i) + g1*( xP2(i) - xP1(i) );
			yPts(IDs(i)) = yP1(i) + g1*( yP2(i) - yP1(i) );	
	
            if (( t > tP1-dt/2) && ( t < tP1+dt/2) )
                xM1 = xPts(IDs);
                yM1 = yPts(IDs);
                please_Print_Vertices_To_File(xP1,yP1,xP2,yP2,xP3,yP3,xM1,yM1);
            end
 
    elseif ((t>tP1)&&(t<=(tP1+tP2)))
			
			%PHASE 2 --> PHASE 2-Better
            
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
            
			xPts(IDs(i)) = xM(i) + g2*( xP2(i) - xM(i) );
			yPts(IDs(i)) = yM(i) + g2*( yP2(i) - yM(i) );            
            
            if (( t > tP1+tP2-dt/2) && ( t < tP1+tP2+dt/2) )
                xM1 = xPts(IDs);
                yM1 = yPts(IDs);
                please_Print_Vertices_To_File(xP1,yP1,xP2,yP2,xP3,yP3,xM1,yM1);
            end
 

            
    elseif ((t>tP1+tP2)&&(t<=(tP1+tP2+tP3)))
			
			%PHASE 2 --> PHASE 2-Even Better!
            
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
            
			xPts(IDs(i)) = xM(i) + g2*( xP2(i) - xM(i) );
			yPts(IDs(i)) = yM(i) + g2*( yP2(i) - yM(i) );
   
    elseif ((t>tP1+tP2+tP3+tP4)&&(t<=(tP1+tP2+tP3+tP4+tP5)))
			
			%PHASE 2 --> PHASE 3!
            
			tprev = tP1+tP2+tP3+tP4;
			t1 = 0.1*tP3 + tprev;
			t2 = 0.9*tP3 + tprev;
			if (t<t1) 							%//For Polynomial Phase Interp.
				g2 = a*power( ( (t-tprev)/tP5) ,2);
            elseif ((t>=t1)&&(t<t2)) 
				g2 = c*power( ( (t-tprev)/tP5) ,3) + d*power( ((t-tprev)/tP5) ,2) + g*( (t-tprev)/tP5) + h;
            elseif (t>=t2) 
				g2 = -b*power( (( (t-tprev)/tP5) - 1) ,2) + 1;
            end			
            
			xPts(IDs(i)) = xP2(i) + g2*( xP3(i) - xP2(i) );
			yPts(IDs(i)) = yP2(i) + g2*( yP3(i) - yP2(i) );
            
%     elseif ((t>tP2)&&(t<=(tP1+tP2+tP3)))
% 			
% 			%PHASE 2 --> PHASE 3
%             
% 			tprev = tP1+tP2;
% 			t1 = 0.1*tP3 + tprev;
% 			t2 = 0.9*tP3 + tprev;
% 			if (t<t1) 							%//For Polynomial Phase Interp.
% 				g2 = a*power( ( (t-tprev)/tP3) ,2);
%             elseif ((t>=t1)&&(t<t2)) 
% 				g2 = c*power( ( (t-tprev)/tP3) ,3) + d*power( ((t-tprev)/tP3) ,2) + g*( (t-tprev)/tP3) + h;
%             elseif (t>=t2) 
% 				g2 = -b*power( (( (t-tprev)/tP3) - 1) ,2) + 1;
%             end			
%             
% 			xPts(IDs(i)) = xP2(i) + g2*( xP3(i) - xP2(i) );
% 			yPts(IDs(i)) = yP2(i) + g2*( yP3(i) - yP2(i) );
%             
    
    end
    
    targets(IDs(i),2) = xPts(IDs(i)); % Store new xVals
    targets(IDs(i),3) = yPts(IDs(i)); % Store new yVals

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints all Vertices to File
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function please_Print_Vertices_To_File(X1,Y1,X2,Y2,X3,Y3,xM,yM)

delete('All_Positions.txt');

size(X1)
size(xM)

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

function [x1,y1,x2,y2,x3,y3,xM,yM] = read_File_In(file_name)

filename = file_name;  %Name of file to read in

fileID = fopen(filename);

    % Read in the file, use 'CollectOutput' to gather all similar data together
    % and 'CommentStyle' to to end and be able to skip lines in file.
    C = textscan(fileID,'%f %f %f %f %f %f %f %f','CollectOutput',1);

fclose(fileID);        %Close the data file.

mat_info = C{1};   %Stores all read in data

%Store all elements in matrix
mat = mat_info(1:end,1:end);

x1 =  mat(:,1); %store regular bingo expectation values 
y1 =  mat(:,2); %store inner bingo expectation values 
x2 =  mat(:,3); %store outer bingo expectation values 
y2 =  mat(:,4); %store cover all bingo expectation values
x3 =  mat(:,5);
y3 =  mat(:,6);
xM =  mat(:,7);
yM =  mat(:,8);