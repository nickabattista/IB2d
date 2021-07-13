%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled  
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  battistn@tcnj.edu
% Date Created: May 27th, 2015
% Institution Created: UNC-CH
% Date Modified: July 13, 2021
% Institution Modified: TCNJ
%
% This code is capable of creating Lagrangian Structures using:
% 	1. Springs
% 	2. Beams (*torsional springs*)
% 	3. Target Points
%	4. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-Tension)")
%
% One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting lengths, etc
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
% 
% If you would like us to add a specific muscle model, please let Nick (battistn@tcnj.edu) know.
%
%--------------------------------------------------------------------------------------------------------------------%


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
