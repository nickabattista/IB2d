%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in the velocity data field from .vtk format
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [U,V] = read_Eulerian_Velocity__Field_vtk(path,simNums)

analysis_path = pwd; % Store path to analysis folder!

cd(path);

filename = ['u.' num2str(simNums) '.vtk'];  % desired lagPts.xxxx.vtk file

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

% Store grid info
N = sscanf(str,'%*s %f %*f %*s',1); 
Nx = N(1);

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
[e_Data,count] = fscanf(fileID,strVec,3*Nx*Nx);
if count ~= 3*Nx*Nx
   error('Problem reading in Eulerian Data.'); 
end

% reshape the matrix into desired data type %
e_Data = reshape(e_Data, 3, count/3); % Reshape (3*Nx*Nx,1) vector to (Nx*Nx,3) matrix
e_Data = e_Data';                     % Store vertices in new matrix

U = e_Data(:,1);       % Store U data
V = e_Data(:,2);       % Store V data

U = reshape(U,Nx,Nx)';  % Reshape (Nx*Nx,1) matrix to  (Nx,Nx)
V = reshape(V,Nx,Nx)';  % Reshape (Nx*Nx,1) matrix to  (Nx,Nx)
 
fclose(fileID);         % Closes the data file.

cd ..;                  % Change directory back to ../viz_IB2d/ directory

cd(analysis_path);      % Change directory back to Data Analysis Folder
