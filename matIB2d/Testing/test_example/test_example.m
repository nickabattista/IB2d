close all
clear 

cd ../Examples

% find all paths to directories and files
list=dir('**/*.*')

% keep paths that lead to directories
counter=0
for i=1:numel(list)
    if(list(i).isdir)>0
        counter=counter+1;
        indx(counter)=(i);
    end
end

fold=list(indx);
counter_id=0;
counter_idv=0;

% iterate over all paths
for i=[3:numel(fold)]
    path=[fold(i).folder '\' fold(i).name]

    % try to cd into path, sometimes the paths don't exist
    try
    % if the path exists cd into the path
    cd(path)
    % check to see if we are in an example directory
    if exist('main2d.m','file')>0
        % try running the main2d.m file if there is an error the path is saved 
        try 
        main2d;
           % if the simulation doesn't create a viz_IB2d folder save the path
           if ~exist('viz_IB2d')
              path
              sprintf('This folder isnt creating viz_IB2d folders')
              counter_idv=counter_idv+1;
              viz_fold(counter_idv)=i;
              fold(i).folder
              continue
           % at the end of the simulation delete the hier_IB2d_data and viz_IB2d directories
           % unless it is one of the restart examples
           elseif (path(end-6:end)=='Restart')
              continue
           else
              if exist('hier_IB2d_data')
                 [status, message, messageid]=rmdir('hier_IB2d_data','s')
              end
              [status, message, messageid] = rmdir('viz_IB2d','s')
        end
        % problem id identifies the examples that give errors
        catch ME
            counter_id=counter_id+1;
            problem_ID(counter_id)=i;
        end
        
        
    end
    catch
    end
end
