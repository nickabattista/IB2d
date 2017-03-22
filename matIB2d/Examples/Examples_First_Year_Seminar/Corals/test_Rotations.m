%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: this function mimics the movement of the simulation for target
%           points based on an angular movement
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function test_Rotations()

angVec = 0:0.05:4*pi;

for i=1:length(angVec)

    % Declare how open you want the movements to go
    full_ang = pi/3;
    dbl_ang = 2*full_ang;

    % Find associated angle
    ang = mod( angVec(i), dbl_ang );
    if ( (ang > full_ang ) && (ang <= dbl_ang) )
        ang = dbl_ang-ang;
    end

    % Read In Pts!
    [xRef,yRef] = read_File_In('hive.vertex');

    plot(xRef,yRef,'b*'); hold on;
    axis([0 1 0 1]);

    % Store Left/Right Points
    xL_Ref = xRef(1:end/2);     yL_Ref = yRef(1:end/2);
    xR_Ref = xRef(end/2+1:end); yR_Ref = yRef(end/2+1:end);

    % Store Values for Centers of Rotation
    xL = xRef(1);       yL = yRef(1);
    xR = xRef(end/2+1); yR = yRef(end/2+1);

    % Rotate Geometry
    [xR_Ref,yR_Ref] = rotate_Geometry(-ang,xR,yR,xR_Ref,yR_Ref);
    [xL_Ref,yL_Ref] = rotate_Geometry(ang,xL,yL,xL_Ref,yL_Ref);

    % Store New Geometry
    xRef2 = [xL_Ref xR_Ref]; 
    yRef2 = [yL_Ref yR_Ref]; 

    plot(xRef2,yRef2,'k*'); hold on;
    axis([0 1 0 1]);
    pause(0.01);
    clf;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: rotate geometry
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,y] = rotate_Geometry(ang,xC,yC,xRef,yRef)

len = length(xRef);
x = zeros(1,len); y=x;

xRef = xRef - xC;
yRef = yRef - yC;

for i=1:len
   x(i) =  xRef(i)*cos(ang) - yRef(i)*sin(ang);
   y(i) =  xRef(i)*sin(ang) + yRef(i)*cos(ang);
end

x = x + xC;
y = y + yC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: reads in info from file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x1,y1] = read_File_In(file_name)

filename = file_name;  %Name of file to read in

fileID = fopen(filename);

    % Read in the file, use 'CollectOutput' to gather all similar data together
    % and 'CommentStyle' to to end and be able to skip lines in file.
    C = textscan(fileID,'%f %f','CollectOutput',1);

fclose(fileID);        %Close the data file.

mat_info = C{1};   %Stores all read in data

%Store all elements in matrix
mat = mat_info(2:end,1:end);

x1 =  mat(:,1); %store xVals 
y1 =  mat(:,2); %store yVals
