%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives (x,y) positions of the immersed boundary at a single step
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag] = give_Lag_Positions(path,numSim)

analysis_path = pwd;

[xLag,yLag] = read_Lagrangian_Data_From_vtk(path,numSim);

cd(analysis_path);

clear analysis_path;
