function test_Script()

% Create Time Vector
tVec=0:0.025:5;

% Loop through time
for i=1:length(tVec) 


% Coefficients for Polynomial Phase-Interpolation
a = 2.739726027397260;  % y1(t) = at^2
b = 2.739726027397260;  % y3(t) = -b(t-1)^2+1
c = -2.029426686960933; % y2(t) = ct^3 + dt^2 + gt + h
d = 3.044140030441400;
g = -0.015220700152207;
h = 0.000253678335870;

% Period Info
tP1 = 0.25;                        % Down-stroke
tP2 = 0.25;                        % Up-stroke
period = tP1+tP2;                   % Period
t = rem(tVec(i),period);       % Current time in simulation ( 'modular arithmetic to get time in period')

% Read In y_Pts for two Phases!
[xP1,yP1,yP2] = read_File_In('swimmer.phases');


%
% FIRST WE COMPUTE THE INTERPOLATE GEOMETRY BETWEEN BOTH PHASES
%

    %PHASE 1 --> PHASE 2
    if (t <= tP1) 

			tprev = 0.0;
			t1 = 0.1*tP1;   
			t2 = 0.9*tP1;
			if (t<t1) 							%For Polynomial Phase Interp.
				g1 = a*power((t/tP1),2);
            elseif ((t>=t1)&&(t<t2)) 
				g1 = c*power((t/tP1),3) + d*power((t/tP1),2) + g*(t/tP1) + h;
            elseif (t>=t2)
				g1 = -b*power(((t/tP1) - 1),2) + 1;
            end
			
            %xPts = xP1 + g1*( xP2 - xP1 );	
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
            
            %xPts = xP2 + g2*( xP1 - xP2 );
			yPts = yP2 + g2*( yP1 - yP2 );
    
    end

   
    plot(xP1,yPts,'.'); hold on;
    axis([0 8 0 8]);
    pause(0.001); 
    clf;
    
end



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

x1 =  mat(:,1);  %store xVals 1/2
y1 =  mat(:,2); %store yVals 1 
y2 =  mat(:,3); %store yVals 2
