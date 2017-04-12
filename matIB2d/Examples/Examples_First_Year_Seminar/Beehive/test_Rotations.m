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
    ang1 = mod( angVec(i), dbl_ang );
    if ( (ang1 > full_ang ) && (ang1 <= dbl_ang) )
        ang1 = dbl_ang-ang1;
    end

    ang2 = mod( angVec(i), dbl_ang );
    if ( (ang2 > full_ang ) && (ang2 <= dbl_ang) )
        ang2 = dbl_ang-ang2;
    end

    % Read In Pts!
    [xRef,yRef] = read_File_In('hive.vertex');

    plot(xRef,yRef,'b*'); hold on;
    axis([0 1 0 1]);

    % Store Left/Right Points
    xL_Ref = xRef(1:end/2);     yL_Ref = yRef(1:end/2);
    xR_Ref = xRef(end/2+1:end); yR_Ref = yRef(end/2+1:end);

    % Store Values for Centers of Rotation
    %xL = xRef(1);       yL = yRef(1);
    %xR = xRef(end/2+1); yR = yRef(end/2+1);

    % Rotate Geometry
    %[xR_Ref,yR_Ref] = rotate_Geometry(-ang,xR,yR,xR_Ref,yR_Ref);
    %[xL_Ref,yL_Ref] = rotate_Geometry(ang,xL,yL,xL_Ref,yL_Ref);
    
    % Store Values for Centers of Rotation
    xL_1 = xRef(1);       yL_1 = yRef(1);           % Left side of Left Pair
    xR_1 = xRef(end/4+1); yR_1 = yRef(end/4+1);     % Right side of Left Pair

    xL_2 = xRef(end/2+1);   yL_2 = yRef(end/2+1);   % Left side of Right Pair
    xR_2 = xRef(3*end/4+1); yR_2 = yRef(3*end/4+1); % Right side of Right Pair

    % -> Rotate Geometry <- %
    % LEFT PAIR %
    [xR_Ref_1,yR_Ref_1] = rotate_Geometry(-ang1,xR_1,yR_1,xRef(end/4+1:end/2),yRef(end/4+1:end/2) );
    [xL_Ref_1,yL_Ref_1] = rotate_Geometry(ang1,xL_1,yL_1,xRef(1:end/4),yRef(1:end/4) );
    % RIGHT PAIR %
    [xR_Ref_2,yR_Ref_2] = rotate_Geometry(-ang2,xR_2,yR_2,xRef(3*end/4+1:end),yRef(3*end/4+1:end) );
    [xL_Ref_2,yL_Ref_2] = rotate_Geometry(ang2,xL_2,yL_2,xRef(end/2+1:3*end/4),yRef(end/2+1:3*end/4) );

    % Store New Geometry
    xRef2 = [xL_Ref_1 xR_Ref_1 xR_Ref_2 xL_Ref_2]; 
    yRef2 = [yL_Ref_1 yR_Ref_1 yR_Ref_2 yL_Ref_2]; 

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
