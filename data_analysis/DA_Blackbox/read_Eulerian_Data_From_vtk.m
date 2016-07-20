%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in the desired Eulerian data from .vtk format
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [e_Data,x,y] = read_Eulerian_Data_From_vtk(path,simNums,strChoice,first)

analysis_path = pwd; % Store path to analysis folder!

cd(path);

filename = [strChoice '.' num2str(simNums) '.vtk'];  % desired lagPts.xxxx.vtk file

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
    y = x';
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
[e_Data,count] = fscanf(fileID,strVec,Nx*Nx);
if count ~= Nx*Nx
   error('\nProblem reading in Eulerian Data.'); 
end

% reshape the matrix into desired data type %
e_Data = reshape(e_Data, Nx, count/Nx); % Reshape vector -> matrix (every 3 entries in vector make into matrix row)
e_Data = e_Data';                       % Store vertices in new matrix

fclose(fileID);                         % Closes the data file.

cd ..;                                  % Change directory back to ../viz_IB2d/ directory

cd(analysis_path);                      % Change directory back to Data Analysis Folder

clear filename fileID str strVec count analysis_path;
