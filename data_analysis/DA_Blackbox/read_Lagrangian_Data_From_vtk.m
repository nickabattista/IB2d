%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in (x,y) positions of the immersed boundary from .vtk
%           format
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag] = read_Lagrangian_Data_From_vtk(path,simNums)

cd(path);

filename = ['lagsPts.' num2str(simNums) '.vtk'];  % desired lagPts.xxxx.vtk file

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

cd ..;                        % Change directory back to ../viz_IB2d/ directory

clear mat str filename fileID count;
