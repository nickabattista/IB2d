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
% FUNCTION: Reads in (x,y) positions of the immersed boundary from .vtk
%           format
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Fdata = read_Force_Scalar_Data_From_vtk(path,simNums,strChoice)

cd(path);

filename = [strChoice '.' num2str(simNums) '.vtk'];  % desired FORCE-LAG-DATA.xxxx.vtk file

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
str = fgets(fileID);

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

for i=1:numLagPts+5 % +5 to get to FORCE DATA
    str = fgets(fileID);
end

% read in the vertices %
[mat,count] = fscanf(fileID,'%f',numLagPts);
if count ~= numLagPts
   error('\nProblem reading in Lagrangian Pts. Force Data'); 
end

mat = reshape(mat, 1, count/1); % Reshape vector -> matrix (every 3 entries in vector make into matrix row)
forces = mat';              % Store vertices in new matrix

fclose(fileID);               % Closes the data file.

Fdata = forces(:,1);         % magnitude of the force

cd ..;                        % Change directory back to ../hier_IB2d_data/ directory

clear mat str filename fileID count;
