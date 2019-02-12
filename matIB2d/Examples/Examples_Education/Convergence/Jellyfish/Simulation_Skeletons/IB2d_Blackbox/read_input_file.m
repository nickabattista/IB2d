%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  nick.battista@unc.edu
% Date Created: May 27th, 2015
% Institution: UNC-CH
%
% This code is capable of creating Lagrangian Structures using:
% 	1. Springs
% 	2. Beams (*torsional springs)
% 	3. Target Points
%   4. Mass Points
%   5. Porous Points
%	6. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-Tension)")
%   7. 3-Element Hill Muscle Model
%
% One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting-lengths, etc
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
% 
% If you would like us to add a specific muscle model, please let Nick (nick.battista@unc.edu) know.
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in the input2d file information
%        
% Author: Aaron Barrett, UNC-CH, abarret@live.unc.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parameters = please_Read_input2d_File(file_name)

file_name = removeComments(file_name);
FID = fopen(file_name, 'r');
parameters = cell(1,2);

tline = fgets(FID);
i = 1;
while(ischar(tline))
    find_brace = strfind(tline,'{');
    if(find_brace)
        parameters{i,1} = strtrim(tline(1:find_brace(1)-1));
        parameters{i,2} = readBrace(FID);
        i = i+1;
    end
    tline = fgets(FID);
end
frewind(FID);
if (isempty(parameters{1,1}))
    % Whoops, no braces
    parameters = readNoBrace(FID);
end
fclose(FID);
delete(file_name);

end

function input_file = removeComments(input_file)
FID = fopen(input_file);
input_file = [input_file, '.temp'];
FID_temp = fopen(input_file, 'w');

tline = fgets(FID);
while(ischar(tline))
    for i = 1:length(tline)
        if(strcmp('%',tline(i)))
            fprintf(FID_temp,'\n');
            break;
        end
        fprintf(FID_temp,'%c',tline(i));
    end
    tline = fgets(FID);
end

fclose(FID);
fclose(FID_temp);
end

function value = readValue(string)
find_quote = strfind(string, '"');
if(find_quote)
    % Value is a string
    value = strtrim(string(find_quote(1)+1:find_quote(2)-1));
else
    % Value is a number
    value = str2num(string);
end
end

function params = readBrace(fid)
params = cell(1,2);
i = 1;
tline = fgets(fid);
brace_close = strfind(tline, '}');
while(isempty(brace_close))
    brace_open = strfind(tline, '{');
    if(brace_open)
        params{i,1} = tline(1:brace_open(1)-1);
        params{i,2} = readBrace(fid);
    else
        find_equal = strfind(tline, '=');
        params{i,1} = tline(1:find_equal(1)-1);
        params{i,2} = readValue(tline(find_equal(1)+1:end));
    end
    tline = fgets(fid);
    brace_close = strfind(tline, '}');
    i = i+1;
end
end

function params = readNoBrace(fid)
params = cell(1,2);
i = 1;
tline = fgets(fid);
while(ischar(tline))
    find_equal = strfind(tline, '=');
    if(find_equal)
        params{i,1} = strtrim(tline(1:find_equal(1)-1));
        params{i,2} = readValue(tline(find_equal(1)+1:end));
        i = i + 1;
    end
    tline = fgets(fid);
end

end