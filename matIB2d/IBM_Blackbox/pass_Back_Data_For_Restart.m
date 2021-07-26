%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  battistn@tcnj.edu
% Date Created: October 8th, 2018
% Institution: TCNJ
%
% This code is capable of creating Lagrangian Structures using:
% 	1. Springs
% 	2. Beams (*torsional springs / non-invariant beams *)
% 	3. Target Points
%	4. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-Tension)")
%
% One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting lengths, etc
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
% 
% If you would like us %to add a specific muscle model, please let Nick (battistn@tcnj.edu) know.
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: provides data from previous .vtk information for restarting a
%           simulation that has ended because of power failure, etc.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [current_time,cter,ctsave,U,V,xLag,yLag,xLag_P,yLag_P] = pass_Back_Data_For_Restart(dt,ctsave,print_dump,path_to_data);


% 
% Automated
%
cter = ctsave * print_dump; % Total # of time-steps thus far up to and included last data point saved
current_time = cter * dt;   % Current time in simulation when last time-step was saved


%
% Read in LAGRANGIAN data for last and second to last timepoint for U,V,mVelocity,F_Poro, xLag,yLag,xLag_P,yLag_P
%
[xLag,yLag,xLag_P,yLag_P] = please_Give_Saved_Lag_Point_Info(ctsave,path_to_data);

%
% Read in EULERIAN data for last timepoint for U,V
%
[U,V] = please_Give_Saved_Eulerian_Info(ctsave,path_to_data);


%
% Update for next iteration (*pretends data from last timepoint was just saved*)
%
ctsave = ctsave + 1;                % Update for next ctsave number (always increases by 1 after data is saved)
current_time = current_time + dt;   % Update for next time-step
cter = cter + 1;                    % Update for next time-step

cd ../ % Go back into Example Directory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: passes back EULERIAN data from last saved timepoint
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [U,V] = please_Give_Saved_Eulerian_Info(ctsave,path_to_data)

    % Points to desired data viz_IB2d data file for CURRENT time-step
    if ctsave<10
       numSim = ['000', num2str(ctsave) ];
    elseif ctsave<100
       numSim = ['00', num2str(ctsave) ];
    elseif ctsave<1000
       numSim = ['0', num2str(ctsave)];
    else
       numSim = num2str(ctsave);
    end
    
    % Imports (x,y) grid values and ALL EULERIAN DATA %
    %                      DEFINITIONS 
    %          x: x-grid                y: y-grid
    %       Omega: vorticity           P: pressure
    %    uMag: mag. of velocity  
    %    uX: mag. of x-Velocity   uY: mag. of y-Velocity  
    %    U: x-directed velocity   V: y-directed velocity
    %    Fx: x-directed Force     Fy: y-directed Force
    %
    %  Note: U(j,i): j-corresponds to y-index, i to the x-index
    %
    % 
    Eulerian_Flags(1) = 0;   % OMEGA
    Eulerian_Flags(2) = 0;   % PRESSURE
    Eulerian_Flags(3) = 0;   % uMAG
    Eulerian_Flags(4) = 0;   % uX (mag. x-component of velocity)
    Eulerian_Flags(5) = 0;   % uY (mag. x-component of velocity)
    Eulerian_Flags(6) = 1;   % uVEC (vector components of velocity: U,V)
    Eulerian_Flags(7) = 0;   % Fx (x-component of force )
    Eulerian_Flags(8) = 0;   % Fy (y-component of force)
    %
    [~,~,~,~,~,~,~,U,V,~,~] = import_Eulerian_Data(path_to_data,numSim,Eulerian_Flags);

    % TURNS OUT UNNECESSARY
    %U = U'; % For IB2d convention from .vtk
    %V = V'; % For IB2d convention from .vtk 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: passes back LAGRANGIAN data from last two saved timepoints
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag,xLag_P,yLag_P] = please_Give_Saved_Lag_Point_Info(ctsave,path_to_data)

    iP = ctsave - 1; % Previous time-point

    % Points to desired data viz_IB2d data file for CURRENT time-step
    if ctsave<10
       numSim = ['000', num2str(ctsave) ];
    elseif ctsave<100
       numSim = ['00', num2str(ctsave) ];
    elseif ctsave<1000
       numSim = ['0', num2str(ctsave)];
    else
       numSim = num2str(ctsave);
    end

    % Points to desired data viz_IB2d data file for PREVIOUS time-step data
    if iP<10
       numSim_Prev = ['000', num2str(iP) ];
    elseif iP<100
       numSim_Prev = ['00', num2str(iP) ];
    elseif iP<1000
       numSim_Prev = ['0', num2str(iP)];
    else
       numSim_Prev = num2str(iP);
    end
    
    % Imports immersed boundary positions %
    [xLag,yLag] = give_Lag_Positions(path_to_data,numSim);
    [xLag_P,yLag_P] = give_Lag_Positions(path_to_data,numSim_Prev);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives (x,y) positions of the immersed boundary at a single step
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag] = give_Lag_Positions(path,numSim)

[xLag,yLag] = read_Lagrangian_Data_From_vtk(path,numSim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in (x,y) positions of the immersed boundary from .vtk
%           format
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag] = read_Lagrangian_Data_From_vtk(path,simNums)


cd(path);

filename = ['lagsPts.' num2str(simNums) '.vtk'];  % desired lagsPts.xxxx.vtk file

fileID = fopen(filename);
if ( fileID== -1 )
    error('\nCANNOT OPEN THE FILE!');
end

str = fgets(fileID); %-1 if eof
if ~strcmp( str(3:5),'vtk');
    error('\nNot in proper VTK format');
end

% read in the header info %
str = fgets(fileID);
str = fgets(fileID);
str = fgets(fileID);
str = fgets(fileID);
str = fgets(fileID); % Up to white space in previous lagsPts.X files (pre June 2021)

% Check whether VTK file has time info. This is a VTK file with time, 
%       need to read the next 3 lines to have read in appropriate
%       number of header lines.
if ~strcmp( str(1), 'P')
	str = fgets(fileID);
    str = fgets(fileID);
    str = fgets(fileID);
end

    
% stores # of Lagrangian Pts. as stated in .vtk file
numLagPts = sscanf(str,'%*s %f %*s',1); 

% read in the vertices %
[mat,count] = fscanf(fileID,'%f %f %f',3*numLagPts);
if count ~= 3*numLagPts
   error('\nProblem reading in Lagrangian Pts.'); 
end

mat = reshape(mat, 3, count/3); % Reshape vector -> matrix (every 3 entries in vector make into matrix row)
vertices = mat';              % Store vertices in new matrix

fclose(fileID);               % Closes the data file.

xLag = vertices(:,1);         % x-Lagrangian Pts.
yLag = vertices(:,2);         % y-Lagrangian Pts.

%---------------------------------------------------------------
% NEEDS TO BE COMMENTED OUT!
% (when compared to other file in data analysis toolbox)
%---------------------------------------------------------------
%cd ..;                        % Change directory back to ../viz_IB2d/ directory

clear mat str filename fileID count;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: imports all Eulerian Data at a single step
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,y,Omega,P,uMag,uX,uY,U,V,Fx,Fy] = import_Eulerian_Data(path,numSim,Eulerian_Flags)

    %
    % EULERIAN FLAGS FOR WHAT GETS SPIT OUT %
    % 
    %Eulerian_Flags(1):   OMEGA
    %Eulerian_Flags(2):   PRESSURE
    %Eulerian_Flags(3):   uMAG
    %Eulerian_Flags(4):   uX (mag. x-component of velocity)
    %Eulerian_Flags(5):   uY (mag. x-component of velocity)
    %Eulerian_Flags(6):   uVEC (x,y-components of velocity: U,V)
    %Eulerian_Flags(7):   Fx (x-component of force )
    %Eulerian_Flags(8):   Fy (y-component of force)
    %
    
%analysis_path = pwd
    
% read in Vorticity %
if Eulerian_Flags(1)
    strChoice = 'Omega'; first = 1;
    [Omega,x,y] = read_Eulerian_Data_From_vtk(path,numSim,strChoice,first);
else
    Omega=[];
end


% read in Pressure %
if Eulerian_Flags(2)
    strChoice = 'P'; first = 1;
    [P,x,y] = read_Eulerian_Data_From_vtk(path,numSim,strChoice,first);
else
    P=[];
end


% read in Velocity Magnitude %
if Eulerian_Flags(3)
    strChoice = 'uMag'; first = 1;
    [uMag,x,y] = read_Eulerian_Data_From_vtk(path,numSim,strChoice,first);
else
    uMag=[];
end

% read in x-directed Velocity Magnitude %
if Eulerian_Flags(4)
    strChoice = 'uX'; first = 1;
    [uX,x,y] = read_Eulerian_Data_From_vtk(path,numSim,strChoice,first);
else
    uX=[];
end

% read in y-directed Velocity Magnitude %
if Eulerian_Flags(5)
    strChoice = 'uY'; first = 1;
    [uY,x,y] = read_Eulerian_Data_From_vtk(path,numSim,strChoice,first);
else
    uY=[];
end

% read in x-directed Forces %
if Eulerian_Flags(7)
    strChoice = 'Fx'; first = 1;
    [Fx,x,y] = read_Eulerian_Data_From_vtk(path,numSim,strChoice,first);
else
    Fx=[];
end

% read in y-directed Forces %
if Eulerian_Flags(8)
    strChoice = 'Fy'; first = 1;
    [Fy,x,y] = read_Eulerian_Data_From_vtk(path,numSim,strChoice,first);
else
    Fy=[];
end

% read in Velocity Field %
if Eulerian_Flags(6)
    [U,V] = read_Eulerian_Velocity_Field_vtk(path,numSim);
else
    U=[];
    V=[];
end

% Default for x,y values
if max(Eulerian_Flags([1:5,7,8]))==0
    x=[];
    y=[];
end

%cd(analysis_path);

clear analysis_path;

clear strChoice first;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in the velocity data field from .vtk format
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s

function [U,V] = read_Eulerian_Velocity_Field_vtk(path,simNums)

analysis_path = pwd; % Store path to analysis folder!

cd(path);

filename = ['u.' num2str(simNums) '.vtk'];  % desired EULERIAN-DATA.xxxx.vtk file

fileID = fopen(filename);
if ( fileID== -1 )
    error('\nCANNOT OPEN THE FILE!');
end

str = fgets(fileID); %-1 if eof
if ~strcmp( str(3:5),'vtk')
    error('\nNot in proper VTK format');
end

% read in the header info %
str = fgets(fileID);
str = fgets(fileID);
str = fgets(fileID);
str = fgets(fileID);
str = fgets(fileID);

% Check whether VTK file has time info. This is a VTK file with time, 
%       need to read the next 3 lines to have read in appropriate
%       number of header lines.
if ~strcmp( str(1:3), 'DIM')
	str = fgets(fileID);
    str = fgets(fileID);
    str = fgets(fileID);
end

% Store grid info
Nx = sscanf(str,'%*s %f %*f %*s',1);
Ny = sscanf(str,'%*s %*f %f %*s',1);


% bypass lines in header %
str = fgets(fileID);
str = fgets(fileID);
str = fgets(fileID);
str = fgets(fileID);
str = fgets(fileID);
str = fgets(fileID);
str = fgets(fileID);


% get formatting for reading in data from .vtk in fscanf %
strVec = '%f';
for i=2:3*Nx
    strVec = [strVec ' %f'];
end

% read in the vertices %
[e_Data,count] = fscanf(fileID,strVec,3*Nx*Ny);
if count ~= 3*Nx*Ny
   error('Problem reading in Eulerian Data.'); 
end

% reshape the matrix into desired data type %
e_Data = reshape(e_Data, 3, count/3); % Reshape (3*Nx*Nx,1) vector to (Nx*Nx,3) matrix
e_Data = e_Data';                     % Store vertices in new matrix

U = e_Data(:,1);       % Store U data
V = e_Data(:,2);       % Store V data

U = reshape(U,Nx,Ny)';  % Reshape (Nx*Nx,1) matrix to  (Nx,Nx)
V = reshape(V,Nx,Ny)';  % Reshape (Nx*Nx,1) matrix to  (Nx,Nx)
 
fclose(fileID);         % Closes the data file.

cd ..;                  % Change directory back to ../viz_IB2d/ directory

cd(analysis_path);      % Change directory back to Data Analysis Folder

clear filename fileID str strVec count analysis_path e_Data Nx;



