function targets = update_Target_Point_Positions(dt,current_time,targets)


IDs = targets(:,1);                 % Stores Lag-Pt IDs in col vector
xPts= targets(:,2);                 % Original x-Values of x-Target Pts.
yPts= targets(:,3);                 % Original y-Values of y-Target Pts.
kStiffs = targets(:,4);             % Stores Target Stiffnesses 

N_target = length(targets(:,1));    %Gives total number of target pts!

for i=1:N_target                    % Loops over all target points!
    
    xPts(IDs(i)) = xPts(IDs(i)) + 0.002;
    
    targets(IDs(i),2) = xPts(IDs(i));
end

