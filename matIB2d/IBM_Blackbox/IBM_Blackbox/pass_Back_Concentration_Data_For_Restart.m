%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  nickabattista@gmail.com
% Date Created: April 26, 2021
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
% If you would like us %to add a specific muscle model, please let Nick (nickabattista@gmail.com) know.
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: provides data from previous .vtk information for restarting a
%           simulation that has ended because of power failure, etc.
%           specifically for Concentration Data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C = pass_Back_Concentration_Data_For_Restart(ctsave,path_to_data)

%
% Read in CONCENTRATION data for last timepoint for Concentration, C
%
C = please_Give_Saved_Concentration_Info(ctsave,path_to_data);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: passes back EULERIAN data from last saved timepoint
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C = please_Give_Saved_Concentration_Info(ctsave,path_to_data)

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
    Eulerian_Flags(6) = 0;   % uVEC (vector components of velocity: U,V)
    Eulerian_Flags(7) = 0;   % Fx (x-component of force )
    Eulerian_Flags(8) = 0;   % Fy (y-component of force)
    Eulerian_Flags(9) = 1;   % Concentration
    %
    [~,~,~,~,~,~,~,~,~,~,~,C] = import_Eulerian_Data(path_to_data,numSim,Eulerian_Flags);

    % TURNS OUT UNNECESSARY
    %U = U'; % For IB2d convention from .vtk
    %V = V'; % For IB2d convention from .vtk 
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: imports all Eulerian Data at a single step
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,y,Omega,P,uMag,uX,uY,U,V,Fx,Fy,C] = import_Eulerian_Data(path,numSim,Eulerian_Flags)

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
    %Eulerian_Flags(9):   Concentration
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

% read in Concentration %
if Eulerian_Flags(9)
    strChoice = 'concentration'; first = 1;
    [C,x,y] = read_Eulerian_Data_From_vtk(path,numSim,strChoice,first);
else
    C=[];
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
% FUNCTION: Reads in the desired Eulerian data from .vtk format
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [e_Data,x,y] = read_Eulerian_Data_From_vtk(path,simNums,strChoice,first)

analysis_path = pwd; % Store path to analysis folder!

cd(path);

filename = [strChoice '.' num2str(simNums) '.vtk'];  % desired EULERIAN-DATA.xxxx.vtk file

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
N = sscanf(str,'%*s %f %f %*s',2);
Nx = N(1); Ny = N(2);

% bypass lines in header %
str = fgets(fileID);
str = fgets(fileID);
str = fgets(fileID);

% Store grid resolution if necessary %
if first == 1
    dx = sscanf(str,'%*s %f %*f %*f',1);
    x=0;
    for i=2:Nx
        x(i) = x(i-1)+dx;
    end
    
    dy = sscanf(str,'%*s %*f %f %*f',1);
    y=0;
    for i=2:Ny
        y(i) = y(i-1)+dy;
    end
    
else
    x=1; y=1; % Store arbitrary values
end

% bypass lines in header %
str = fgets(fileID);
str = fgets(fileID);
str = fgets(fileID);
str = fgets(fileID);
str = fgets(fileID);

% get formatting for reading in data from .vtk in fscanf %
strVec = '%f';
for i=2:Nx
    strVec = [strVec ' %f'];
end

% read in the vertices %
[e_Data,count] = fscanf(fileID,strVec,Nx*Ny);
if count ~= Nx*Ny
   error('\nProblem reading in Eulerian Data.'); 
end

% reshape the matrix into desired data type %
e_Data = reshape(e_Data, Nx, count/Nx); % Reshape vector -> matrix (every 3 entries in vector make into matrix row)
e_Data = e_Data';                       % Store vertices in new matrix

fclose(fileID);                         % Closes the data file.

cd ..;                                  % Change directory back to ../viz_IB2d/ directory

cd(analysis_path);                      % Change directory back to Data Analysis Folder

clear filename fileID str strVec count analysis_path;
