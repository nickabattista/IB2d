%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  nickabattista@gmail.com
% Date Created: May 27th, 2015
% Institution: UNC-CH
%
% This code is capable of creating Lagrangian Structures using:
% 	1. Springs
% 	2. Beams (torsional springs and non-invariant beams)
% 	3. Target Points
%	4. Muscle-Model (combined Force-Length-Velocity model, "Hill+(Length-Tension)")
%
% One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting lengths, etc
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
% 
% If you would like us to add a specific muscle model, please let Nick (nickabattista@gmail.com) know.
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: FUNCTION computes the Lagrangian COAGULATION Force Densities
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fx, fy, aggregate_list] = give_Me_Coagulation_Force_Densities(Nb,xLag,yLag,coagulation,aggregate_list,Lx,Ly)


    % 1st: Find Cell (x,y) Centers
    % 2nd: Find Distances Between Cells and Make List of Close Connections
    % 3rd: Find Distances of What Lag. Pts. Are Close Enough to Near Center
    % 4th: Make Spring Connections (resting_length = original nearest distance), e.g., the bond
    % 5th: Check Spring Connections for Threshold Forces and Remove Spring Connection, e.g., remove the bond
    
    % coagulation:     row 1: staring index for coagulation cells
    %                  row 2: threshold radii
    %                  row 3: bond strength
    %                  row 4: threshold force
    %                  row 5: row N+1: # of Lag. Pts in Cell, 1
    %                  row 6: row N+1: # of Lag. Pts in Cell, 2
    %                   .  
    %                   .  
    %                   . 
    %                  row N+4: # of Lag. Pts in Cell, N
    
    start_id = coagulation(1);
    threshold_radius = coagulation(2);
    bond_strength = coagulation(3);
    fracture_force = coagulation(4);
    
    % Find Centers of each Cell
    [xC, yC, Ncells] = please_Find_Cell_Centers(coagulation(5:end),start_id,xLag,yLag);
    
    % Compute distance between each cell and store connection if less than threshold AND NOT already a bond between cells
    coag_connects = please_Find_Distances_Between_Cells(xC,yC,Ncells,threshold_radius,aggregate_list);

    % find indices for placing NEW bond between
    if coag_connects > 0
        % Find closest lag-pts for each cell
        lag_indices_connects = please_Find_Closest_Lag_Pts_To_Cell_Centers(start_id,xC,yC,xLag,yLag,coag_connects,coagulation(5:end),threshold_radius);
    end

    
    % Make Spring Connections and Compute Forces!
    if coag_connects > 0
        if (aggregate_list ~= 0)
            Ncoag_already = length( aggregate_list(:,1) );                           % # of bonds already previous to this time-step
            aggregate_list = [aggregate_list; coag_connects lag_indices_connects];   % Store all info for connections (NOTE: RL_Vec previously incoroporated into lag_indices_connects matrix)
            [fx, fy, remove_Bonds] = give_Me_Coagulation_Spring_Lagrangian_Force_Densities(Nb,xLag,yLag,aggregate_list,Ncoag_already,fracture_force,bond_strength,Lx,Ly);
        else
            Ncoag_already = 0;
            aggregate_list = [coag_connects lag_indices_connects];   % Store all info for connections (NOTE: RL_Vec previously incoroporated into lag_indices_connects matrix)
            [fx, fy, remove_Bonds] = give_Me_Coagulation_Spring_Lagrangian_Force_Densities(Nb,xLag,yLag,aggregate_list,Ncoag_already,fracture_force,bond_strength,Lx,Ly);
        end
    else
         if aggregate_list ~= 0
            Ncoag_already = length( aggregate_list(:,1) );                          % # of bonds already previous to this time-step
            [fx, fy, remove_Bonds] = give_Me_Coagulation_Spring_Lagrangian_Force_Densities(Nb,xLag,yLag,aggregate_list,Ncoag_already,fracture_force,bond_strength,Lx,Ly);
         else
             remove_Bonds = 0;
             fx = zeros(Nb,1);
             fy = zeros(Nb,1);
         end
    end

    
    % Remove Bonds from aggregate_list (based off row marker indices in remove_Bonds)
    if remove_Bonds ~= 0
        aggregate_list = please_Remove_Bonds(remove_Bonds,aggregate_list);
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes center (x,y)-coordinate for each cell
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

function [xC, yC, Ncells] = please_Find_Cell_Centers(coagulation,start_id,xLag,yLag)

    Ncells = length(coagulation);   % # of Cells

    xC = zeros(Ncells,1);
    yC = xC;
    
    for i=1:Ncells
       
        num_Pts_Cell = coagulation(i);  % # of Lag Pts. on Cell, i
        
        % Gives index of first Lag. Pt.
        if i>1
            lag1_ind = start_id + sum(coagulation(1:i-1));
        else
            lag1_ind = start_id;
        end
        
        % Gives index of 2nd Lag. Pt.
        lag2_ind = lag1_ind + 0.5*num_Pts_Cell; 
        
        % Gives Center
        xC(i) = 0.5*( xLag(lag1_ind) + xLag(lag2_ind) );
        yC(i) = 0.5*( yLag(lag1_ind) + yLag(lag2_ind) );
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes distances between cells
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

function connects = please_Find_Distances_Between_Cells(xC,yC,Ncells,threshold_radius,aggregate_list)

connects = 0;   % matrix of NEW connections
num = 1;        % initialize counter
first = 1;      % initalize for storing aggregate_list(:,2)

for i=1:Ncells         % Loops over every cell
   for j=i+1:Ncells      % Loops over all cells not tested against it yet
      
       % Compute distance between every cell
       r = sqrt( ( xC(i) - xC(j) )^2 + ( yC(i) - yC(j) )^2 );
      
       % Check is r is less than threshold
       if r < threshold_radius 
           
           % Check if no current bonds between cells
           if aggregate_list == 0 
               connects(num,1) = i;
               connects(num,2) = j;
               num=num+1;
               
           % If there are already some bonds, check to make sure bond between i-j doesn't already exist, if not, add one.
           else
               
               % Store aggregate 2nd column if first pass thru
               if first == 1
                   agg_list2 = aggregate_list(:,2);
                   first = 0;
               end
               
               % Finds indices associated with i in aggregate_list already
               inds_i = find( aggregate_list(:,1)== i);
               
               % if no indices match i in aggregate_list, store i-j connection
               if inds_i == 0
                   connects(num,1) = i;
                   connects(num,2) = j;
                   num=num+1;
               
               % if indices match i in aggregate list 1st column, check if already a connection
               else
                   checkVec = agg_list2(inds_i);
                   if find(checkVec == j)
                        % do nothing, already a bond between them
                   else
                      connects(num,1) = i;
                      connects(num,2) = j;
                      num=num+1;
                   end
               end
               
           end
       
       end
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: finds closest Lagrangian indices to one another
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function lag_indices_connects = please_Find_Closest_Lag_Pts_To_Cell_Centers(start_id,xC,yC,xLag,yLag,connects,coagulation,threshold_radius)

Nconnects = length(connects(:,1));          % # of new bonds 
lag_indices_connects = zeros(Nconnects,2);  % initialize storage for indices
RL_Vec = zeros(Nconnects,1);                % initialize storage for resting-lengths

for i=1:Nconnects
    
    C_ind_1 = connects(i,1); % index of cell 1
    C_ind_2 = connects(i,2); % index of cell 2
    
    st_Cell1_Lag_inds = start_id + sum( coagulation(1:C_ind_1-1) ); % starting Lagrangian index on cell
    st_Cell2_Lag_inds = start_id + sum( coagulation(1:C_ind_2-1) ); % starting Lagrangian index on cell

    dist = threshold_radius+1; % initialize larger than threshold radius
    for j=0:coagulation(C_ind_1)-1
        
        % index associated with Lag. Pt. along Cell, C_ind_1
        ind = st_Cell1_Lag_inds + j;
        
        % dist between that Lag. Pt. and Center of other Cell
        chk = sqrt( ( xLag(ind) - xC(C_ind_2) )^2 + ( yLag(ind) - yC(C_ind_2))^2  );
        
        % check if minimal so far and store index of Lag. Pt.
        if chk < dist
            dist = chk;
            lag_indices_connects(i,1) = ind;
        end
        
    end
    
    dist = threshold_radius+1;   % initialize larger than threshold radius
    for j=0:coagulation(C_ind_2)-1
        
        % index associated with Lag. Pt. along Cell, C_ind_2
        ind = st_Cell2_Lag_inds + j;
        
        % dist between that Lag. Pt. and Center of other Cell
        chk = sqrt( ( xLag(ind) - xC(C_ind_1) )^2 + ( yLag(ind) - yC(C_ind_1))^2  );
        
        % check if minimal so far and store index of Lag. Pt.
        if chk < dist
            dist = chk;
            lag_indices_connects(i,2) = ind;
        end
        
    end
    
    RL_Vec(i,1) = sqrt( ( xLag(lag_indices_connects(i,1)) - xLag(lag_indices_connects(i,2)) )^2 + ( yLag(lag_indices_connects(i,1)) - yLag(lag_indices_connects(i,2)) )^2  );
       
end % End loop over all new bonds


% Combine into one matrix
lag_indices_connects = [lag_indices_connects RL_Vec];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION computes the Lagrangian SPRING Force Densities.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fx, fy, remove_Bonds] = give_Me_Coagulation_Spring_Lagrangian_Force_Densities(Nb,xLag,yLag,agg_list,Ncoag_already,fracture_force,bond_strength,Lx,Ly)


Nsprings = length(agg_list(:,1)); % # of Springs (BONDS)
sp_1 = agg_list(:,3);             % Initialize storage for MASTER NODE Spring Connection
sp_2 = agg_list(:,4);             % Initialize storage for SLAVE NODE Spring Connection
RL_Vec = agg_list(:,5);           % Resting length of bond (whatever initial distance is when bond is formed)
alpha_pow = 1;                    % Degree of linearity (1=linear, >1 = non-linear)

fx = zeros(Nb,1);                 % Initialize storage for x-forces
fy = fx;                          % Initialize storage for y-forces

num = 1;                          % Initialize counter for bond removal
remove_Bonds = 0;                 % Initialize for return

for i=1:Nsprings
    
    id_Master = sp_1(i);          % Master Node index
    id_Slave = sp_2(i);           % Slave Node index
    k_Spring = bond_strength;     % Spring stiffness of i-th spring (CONSTANT FOR NOW)
    L_r = RL_Vec(i);              % Resting length of i-th spring
    alpha = alpha_pow;            % Degree of linearity of i-th spring (CONSTANT FOR NOW)
 
    dx = xLag(id_Slave) - xLag(id_Master); % x-Distance btwn slave and master node
    dy = yLag(id_Slave) - yLag(id_Master); % y-Distance btwn slave and master node
    
    %
    % TESTING FOR LAG PT. PASSED THRU BNDRY; MAY NEED TO CHANGE TOLERANCE HERE, DEPENDENT ON APPLICATION
    %
    if abs(dx) > Lx/2
        dx = sign(dx)*( Lx - sign(dx)*dx );
    end
    
    if abs(dy) > Ly/2
        dy = sign(dy)*( Ly - sign(dy)*dy );
    end
    
    % Compute spring deformation force
    sF_x = 0.5*(alpha+1) * k_Spring * ( sqrt( dx^2 + dy^2 ) - L_r )^(alpha) * ( dx / sqrt(dx^2+dy^2) );
    sF_y = 0.5*(alpha+1) * k_Spring * ( sqrt( dx^2 + dy^2 ) - L_r )^(alpha) * ( dy / sqrt(dx^2+dy^2) );
    
    % CHECK IF PREVIOUS BONDS SHOULD BE BROKEN
    if i <= Ncoag_already
        
        fMag = ( sF_x^2 + sF_y^2 );
        
        if fMag < fracture_force

            fx(id_Master,1) = fx(id_Master,1) + sF_x;  % Sum total forces for node, i in x-direction (this is MASTER node for this spring)
            fy(id_Master,1) = fy(id_Master,1) + sF_y;  % Sum total forces for node, i in y-direction (this is MASTER node for this spring)

            fx(id_Slave,1) = fx(id_Slave,1) - sF_x;    % Sum total forces for node, i in x-direction (this is SLAVE node for this spring)
            fy(id_Slave,1) = fy(id_Slave,1) - sF_y;    % Sum total forces for node, i in y-direction (this is SLAVE node for this spring)

        else

            remove_Bonds(num) = i; % store what row to remove from agg_list
            num=num+1;             % iterate
            
        end
    
    % STORE FORCE FOR NEW BOND
    else
        
        fx(id_Master,1) = fx(id_Master,1) + sF_x;  % Sum total forces for node, i in x-direction (this is MASTER node for this spring)
        fy(id_Master,1) = fy(id_Master,1) + sF_y;  % Sum total forces for node, i in y-direction (this is MASTER node for this spring)

        fx(id_Slave,1) = fx(id_Slave,1) - sF_x;    % Sum total forces for node, i in x-direction (this is SLAVE node for this spring)
        fy(id_Slave,1) = fy(id_Slave,1) - sF_y;    % Sum total forces for node, i in y-direction (this is SLAVE node for this spring)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION removes bonds where force has fractured the bond between cells
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function aggregate_list = please_Remove_Bonds(remove_Bonds,aggregate_list)

%aggregate_list;

% place zeros as markers to get rid of spots
aggregate_list(remove_Bonds,:) = 0;

% only non-zero entries
aggregate_list = aggregate_list(find(aggregate_list > 0));

% reshape into previous form 
%aggregate_list;
N = length(aggregate_list);   % # of total entries in vector formed by find above
aggregate_list = reshape(aggregate_list,N/5,5);