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
% 	2. Beams (*torsional springs)
% 	3. Target Points
%	4. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-Tension)")
%   5. etc, etc.
%
% One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting lengths, etc
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
% 
% If you would like us to add a specific muscle model, please let Nick (nickabattista@gmail.edu) know.
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: reads in input2d files and initializes the simulation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Fluid_Params, Grid_Params, Time_Params, Lag_Struct_Params, Output_Params, Lag_Name_Params,Con_Params] = please_Initialize_Simulation()

% Please initialize all parameters to off (as a cautionary check) %
% [Lag_Struct_Params,Output_Params] = please_Initialize_Off();

% READ IN ALL INPUTS INTO CELLS FROM INPUT2D %
params = please_Read_input2d_File('input2d');

% EXTRACT INDIVIDUAL CELL GROUPS %
Fluid_Input = params{find(strcmp({params{:,1}},'Fluid_Parameters')),2};
Grid_Input = params{find(strcmp({params{:,1}},'Grid_Parameters')),2};
Time_Input = params{find(strcmp({params{:,1}},'Temporal_Information')),2};
Lag_Struct_Input = params{find(strcmp({params{:,1}},'Lag_Structure_Info')),2};
Output_Input = params{find(strcmp({params{:,1}},'Output_Info')),2};
Lag_Name_Input = params{find(strcmp({params{:,1}},'Lag_Name')),2};

try
Con_Input = params{find(strcmp({params{:,1}},'Concentration_Info')),2};
catch
'No Concentration'
Con_Input=cell(2,6);
end
% INITIALIZE PARAMETERS FOR IBM_DRIVER FILE %
Fluid_Params = please_Initialize_Fluid_Inputs(Fluid_Input);
Grid_Params = please_Initialize_Grid_Inputs(Grid_Input);
Time_Params = please_Initialize_Time_Inputs(Time_Input);
Lag_Struct_Params = please_Initialize_Lag_Structure_Inputs(Lag_Struct_Input);
Con_Params = please_Initialize_Con_Inputs(Con_Input);
Output_Params = please_Initialize_Output_Inputs(Output_Input);
Lag_Name_Params = please_Initialize_Lag_Name_Inputs(Lag_Name_Input);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: initializes LAGRANGIAN NAME parameters for IBM_Driver file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Lag_Name_Params = please_Initialize_Lag_Name_Inputs(Lag_Name_Input)

% Lag_Name_Params: string_name

try
    Lag_Name_Params = Lag_Name_Input{find(strcmp({Lag_Name_Input{:,1}},'string_name')),2};
catch
  fprintf('\n \n');
  error('LAGRANGIAN NAME Parameter Improperly Declared in input2d file');   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: initializes FLUID parameters for IBM_Driver file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Fluid_Params = please_Initialize_Fluid_Inputs(Fluid_Input)

% Fluid_Params(1): mu
%             (2): density

try
    Fluid_Params(1) = Fluid_Input{find(strcmp({Fluid_Input{:,1}},'mu')),2};
    Fluid_Params(2) = Fluid_Input{find(strcmp({Fluid_Input{:,1}},'rho')),2};
catch
  fprintf('\n \n');
  error('FLUID Parameters Improperly Declared in input2d file');   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: initializes GRID parameters for IBM_Driver file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Grid_Params = please_Initialize_Grid_Inputs(Grid_Input)

% Grid_Params(1): Nx
%            (2): Ny
%            (3): Lx
%            (4): Ly
%            (5): Supp

try
    Grid_Params(1) = Grid_Input{find(strcmp({Grid_Input{:,1}},'Nx')),2};
    Grid_Params(2) = Grid_Input{find(strcmp({Grid_Input{:,1}},'Ny')),2};
    Grid_Params(3) = Grid_Input{find(strcmp({Grid_Input{:,1}},'Lx')),2};
    Grid_Params(4) = Grid_Input{find(strcmp({Grid_Input{:,1}},'Ly')),2};
    Grid_Params(5) = Grid_Input{find(strcmp({Grid_Input{:,1}},'supp')),2};
catch
    fprintf('\n \n');
    error('GRID Information Improperly Declared in input2d file'); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: initializes TIME parameters for IBM_Driver file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Time_Params = please_Initialize_Time_Inputs(Time_Input)

% Time_Params(1): Tfinal (end time of simulation)
%            (2): dt (time-step)

try
    Time_Params(1) = Time_Input{find(strcmp({Time_Input{:,1}},'Tfinal')),2};
    Time_Params(2) = Time_Input{find(strcmp({Time_Input{:,1}},'dt')),2};
    
    % RESTART FLAG %
    if find(strcmp({Time_Input{:,1}},'Restart_Flag')) > 0
        Time_Params(3) = Time_Input{find(strcmp({Time_Input{:,1}},'Restart_Flag')),2};
    else
        Time_Params(3) = 0;
    end
    
catch
   fprintf('\n \n');
   error('TEMPORAL Information Improperly Declared in input2d file'); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: initializes LAGRANGIAN STRUCTURE parameters for IBM_Driver file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Lag_Struct_Params = please_Initialize_Lag_Structure_Inputs(Lag_Struct_Input)

% Lag_Struct_Params(1): springs
%                  (2): update_springs
%                  (3): target points
%                  (4): update_target_points
%                  (5): beams (torsional beams)
%                  (6): update_beams
%                   .         .
%                   .         .
%                   .         .

try

    % SPRINGS %
    if find(strcmp({Lag_Struct_Input{:,1}},'springs')) > 0
        Lag_Struct_Params(1) = Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'springs')),2};
    else
        Lag_Struct_Params(1) = 0;
    end

    % UPDATE_SPRINGS %
    if find(strcmp({Lag_Struct_Input{:,1}},'update_springs')) > 0
        Lag_Struct_Params(2) = Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'update_springs')),2};
    else
        Lag_Struct_Params(2) = 0;
    end

    % TARGET POINTS %
    if find(strcmp({Lag_Struct_Input{:,1}},'target_pts')) > 0
        Lag_Struct_Params(3) = Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'target_pts')),2};
    else
        Lag_Struct_Params(3) = 0;
    end

    % UPDATE TARGET POINTS %
    if find(strcmp({Lag_Struct_Input{:,1}},'update_target')) > 0
        Lag_Struct_Params(4) = Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'update_target')),2};
    else
        Lag_Struct_Params(4) = 0;
    end

    % BEAMS (torsional springs) %
    if find(strcmp({Lag_Struct_Input{:,1}},'beams')) > 0
        Lag_Struct_Params(5) = Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'beams')),2};
    else
        Lag_Struct_Params(5) = 0;
    end

    % UPDATE BEAMS (torsional springs) %
    if find(strcmp({Lag_Struct_Input{:,1}},'update_beams')) > 0
        Lag_Struct_Params(6) = Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'update_beams')),2};
    else
        Lag_Struct_Params(6) = 0;
    end


    % NON-INVARIANT BEAMS %
    if find(strcmp({Lag_Struct_Input{:,1}},'nonInvariant_beams')) > 0
        Lag_Struct_Params(7) = Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'nonInvariant_beams')),2};
    else
        Lag_Struct_Params(7) = 0;
    end

    % UPDATE NON-INVARIANT BEAMS %
    if find(strcmp({Lag_Struct_Input{:,1}},'update_nonInv_beams')) > 0
        Lag_Struct_Params(8) = Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'update_nonInv_beams')),2};
    else
        Lag_Struct_Params(8) = 0;
    end


    % FV-LT MUSCLE MODEL %
    if find(strcmp({Lag_Struct_Input{:,1}},'FV_LT_muscle')) > 0
        Lag_Struct_Params(9) = Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'FV_LT_muscle')),2};
    else
        Lag_Struct_Params(9) = 0;
    end

    % 3-ELEMENT HILL MODEL %
    if find(strcmp({Lag_Struct_Input{:,1}},'3_element_muscle')) > 0
        Lag_Struct_Params(10) = Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'3_element_muscle')),2};
    else
        Lag_Struct_Params(10) = 0;
    end

    % ARBITRARY EXTERNAL FORCE %
    if find(strcmp({Lag_Struct_Input{:,1}},'arb_ext_force')) > 0
        Lag_Struct_Params(11) = Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'arb_ext_force')),2};
    else
        Lag_Struct_Params(11) = 0;
    end

    % TRACERS %
    if find(strcmp({Lag_Struct_Input{:,1}},'tracers')) > 0
        Lag_Struct_Params(12) = Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'tracers')),2};
    else
        Lag_Struct_Params(12) = 0;
    end

    % MASSIVE POINTS %
    if find(strcmp({Lag_Struct_Input{:,1}},'mass_pts')) > 0
        Lag_Struct_Params(13) = Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'mass_pts')),2};
    else
        Lag_Struct_Params(13) = 0;
    end

    % GRAVITY %
    if find(strcmp({Lag_Struct_Input{:,1}},'gravity')) > 0
        Lag_Struct_Params(14) = Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'gravity')),2};
    else
        Lag_Struct_Params(14) = 0;
    end

    % x-GRAVITY VECTOR COMPONENT %
    if ( find(strcmp({Lag_Struct_Input{:,1}},'x_gravity_vec_comp')) > 0 || find(strcmp({Lag_Struct_Input{:,1}},'x_gravity_vec_comp')) < 0 )
        Lag_Struct_Params(15) = Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'x_gravity_vec_comp')),2};
    else
        Lag_Struct_Params(15) = 0;
    end

    % y-GRAVITY VECTOR COMPONENT %
    if ( find(strcmp({Lag_Struct_Input{:,1}},'y_gravity_vec_comp')) > 0 || find(strcmp({Lag_Struct_Input{:,1}},'y_gravity_vec_comp')) < 0 )
        Lag_Struct_Params(16) = Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'y_gravity_vec_comp')),2};
    else
        Lag_Struct_Params(16) = 0;
    end

    % POROUS MEDIA %
    if find(strcmp({Lag_Struct_Input{:,1}},'porous_media')) > 0
        Lag_Struct_Params(17) = Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'porous_media')),2};
    else
        Lag_Struct_Params(17) = 0;
    end

    % ELECTROPHYSIOLOGY %
    if find(strcmp({Lag_Struct_Input{:,1}},'electro_phys')) > 0
        Lag_Struct_Params(18) = Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'electro_phys')),2};
    else
        Lag_Struct_Params(18) = 0;
    end

    % DAMPED SPRINGS %
    if find(strcmp({Lag_Struct_Input{:,1}},'damped_springs')) > 0
        Lag_Struct_Params(19) = Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'damped_springs')),2};
    else
        Lag_Struct_Params(19) = 0;
    end


    % UPDATE DAMPED SPRINGS %
    if find(strcmp({Lag_Struct_Input{:,1}},'update_damp_springs')) > 0
        Lag_Struct_Params(20) = Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'update_damp_springs')),2};
    else
        Lag_Struct_Params(20) = 0;
    end

    % BOUSSINESQ %
    if find(strcmp({Lag_Struct_Input{:,1}},'boussinesq')) > 0
        Lag_Struct_Params(21) = Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'boussinesq')),2};
    else
        Lag_Struct_Params(21) = 0;
    end

    % EXPANSION COEFFICIENT (BOUSSINESQ) %
    if find(strcmp({Lag_Struct_Input{:,1}},'expansion_coeff')) > 0
        Lag_Struct_Params(22) = Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'expansion_coeff')),2};
    else
        Lag_Struct_Params(22) = 0;
    end

    % USER-FORCE MODEL %
    if find(strcmp({Lag_Struct_Input{:,1}},'user_force_model')) > 0
        Lag_Struct_Params(23) = Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'user_force_model')),2};
    else
        Lag_Struct_Params(23) = 0;
    end
    
    % PORO-ELASTIC MEDIA %
    if find(strcmp({Lag_Struct_Input{:,1}},'poroelastic')) > 0
        Lag_Struct_Params(24) = Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'poroelastic')),2};
    else
        Lag_Struct_Params(24) = 0;
    end
    
    % COAGULATION %
    if find(strcmp({Lag_Struct_Input{:,1}},'coagulation')) > 0
        Lag_Struct_Params(25) = Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'coagulation')),2};
    else
        Lag_Struct_Params(25) = 0;
    end
    
    % BRINKMAN %
    if find(strcmp({Lag_Struct_Input{:,1}},'brinkman')) > 0
        Lag_Struct_Params(26) = Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'brinkman')),2};
    else
        Lag_Struct_Params(26) = 0;
    end

catch
    fprintf('\n \n');
    error('LAGRANGIAN STRUCTURE INFO Improperly Declared in input2d file');     
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: initializes CONCENTRATION parameters for IBM_Driver file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Con_Params = please_Initialize_Con_Inputs(Con_Input)

% Con_Params(1): concentration
%           (2): kDiff
%           (3): advection
%           (4): source
%           (5): k_source
%           (6): c_inf

try

    % CONCENTRATION GRADIENT %
    if find(strcmp({Con_Input{:,1}},'concentration')) > 0
        Con_Params(1) = Con_Input{find(strcmp({Con_Input{:,1}},'concentration')),2};
    else
        Lag_Struct_Params(1) = 0;
    end

    % DIFFUSION COEFFICIENT %
    if find(strcmp({Con_Input{:,1}},'kDiff')) > 0
        Con_Params(2) = Con_Input{find(strcmp({Con_Input{:,1}},'kDiff')),2};
    else
        Con_Params(2) = 0;
    end
    % CONCENTRATION GRADIENT %
    if find(strcmp({Con_Input{:,1}},'advection')) > 0
        Con_Params(3) = Con_Input{find(strcmp({Con_Input{:,1}},'advection')),2};
    else
        Con_Params(3) = 0;
    end


    % CONCENTRATION SOURCE TERM %
    if find(strcmp({Con_Input{:,1}},'source')) > 0
        Con_Params(4) = Con_Input{find(strcmp({Con_Input{:,1}},'source')),2};
    else
        Con_Params(4) = 0;
    end

    % CONCENTRATION SOURCE TERM %
    if find(strcmp({Con_Input{:,1}},'k_source')) > 0
        Con_Params(5) = Con_Input{find(strcmp({Con_Input{:,1}},'k_source')),2};
    else
        Con_Params(5) = 0;
    end

    % CONCENTRATION SATURATION LIMIT %
    if find(strcmp({Con_Input{:,1}},'c_inf')) > 0
        Con_Params(6) = Con_Input{find(strcmp({Con_Input{:,1}},'c_inf')),2};
    else
        Con_Params(6) = 0;
    end
    
    % CONCENTRATION PERIODIC BOUNDARY % 
    %if find(strcmp({Con_Input{:,1}},'periodic_boundary')) > 0
    %    Con_Params(7) = Con_Input{find(strcmp({Con_Input{:,1}},'periodic_boundary')),2};
    %else
    %    Con_Params(7) = 0;
    %end

catch
    fprintf('\n \n');
    error('CONCENTRATION INFO Improperly Declared in input2d file');     
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: initializes OUTPUT parameters for IBM_Driver file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Output_Params = please_Initialize_Output_Inputs(Output_Input)

% Output_Params(1):  print_dump
%              (2):  plot_Matlab
%              (3):  plot_LagPts
%              (4):  plot_Velocity
%              (5):  plot_Vorticity
%              (6):  plot_MagVelocity
%              (7):  plot_Pressure
%              (8):  plot_Concentration
%              (9):  save_Vorticity 
%              (10):  save_Pressure 
%              (11): save_uVec 
%              (12): save_uMag 
%              (13): save_uX 
%              (14): save_uY 
%              (15): save_fMag 
%              (16): save_fX 
%              (17): save_fY 
%              (18): save_hier 

try
    
    % PRINT DUMP INTERVAL %
    if find(strcmp({Output_Input{:,1}},'print_dump')) > 0
        Output_Params(1) = Output_Input{find(strcmp({Output_Input{:,1}},'print_dump')),2};
    elseif find(strcmp({Output_Input{:,1}},'print_Dump')) > 0
        Output_Params(1) = Output_Input{find(strcmp({Output_Input{:,1}},'print_Dump')),2};
    else
        Output_Params(1) = 0;
    end

    %%%%%%                                         %%%%%%
    %%%%%% OUTPUT INFO PLOTTING DIRECTLY IN MATLAB %%%%%%
    %%%%%%                                         %%%%%%
    
    % PLOT IN MATLAB FLAG %
    if find(strcmp({Output_Input{:,1}},'plot_Matlab')) > 0
        Output_Params(2) = Output_Input{find(strcmp({Output_Input{:,1}},'plot_Matlab')),2};
    else
        Output_Params(2) = 0;
    end

    % PLOT LAG. PTS IN MATLAB %
    if find(strcmp({Output_Input{:,1}},'plot_LagPts')) > 0
        Output_Params(3) = Output_Input{find(strcmp({Output_Input{:,1}},'plot_LagPts')),2};
    else
        Output_Params(3) = 0;
    end

    % PLOT VELOCITY IN MATLAB %
    if find(strcmp({Output_Input{:,1}},'plot_Velocity')) > 0
        Output_Params(4) = Output_Input{find(strcmp({Output_Input{:,1}},'plot_Velocity')),2};
    else
        Output_Params(4) = 0;
    end

    % PLOT VORTICITY IN MATLAB %
    if find(strcmp({Output_Input{:,1}},'plot_Vorticity')) > 0
        Output_Params(5) = Output_Input{find(strcmp({Output_Input{:,1}},'plot_Vorticity')),2};
    else
        Output_Params(5) = 0;
    end


    % PLOT MAG. VELOCITY IN MATLAB %
    if find(strcmp({Output_Input{:,1}},'plot_MagVelocity')) > 0
        Output_Params(6) = Output_Input{find(strcmp({Output_Input{:,1}},'plot_MagVelocity')),2};
    else
        Output_Params(6) = 0;
    end
    
    % PLOT PRESSURE IN MATLAB %
    if find(strcmp({Output_Input{:,1}},'plot_Pressure')) > 0
        Output_Params(7) = Output_Input{find(strcmp({Output_Input{:,1}},'plot_Pressure')),2};
    else
        Output_Params(7) = 0;
    end

    % PLOT CONCENTRATION IN MATLAB %
    if find(strcmp({Output_Input{:,1}},'plot_Concentration')) > 0
        Output_Params(8) = Output_Input{find(strcmp({Output_Input{:,1}},'plot_Concentration')),2};
    else
        Output_Params(8) = 0;
    end
    
    %%%%%%                                                %%%%%%
    %%%%%% OUTPUT INFO FOR WHAT GETS SAVED TO .VTK FORMAT %%%%%%
    %%%%%%                                                %%%%%%
    
    % SAVE VORTICITY TO .VTK FORMAT %
    if find(strcmp({Output_Input{:,1}},'save_Vorticity')) >= 0
        Output_Params(9) = Output_Input{find(strcmp({Output_Input{:,1}},'save_Vorticity')),2};
    else
        Output_Params(9) = 1;
    end
    
    % SAVE PRESSURE TO .VTK FORMAT %
    if find(strcmp({Output_Input{:,1}},'save_Pressure')) >= 0
        Output_Params(10) = Output_Input{find(strcmp({Output_Input{:,1}},'save_Pressure')),2};
    else
        Output_Params(10) = 1;
    end
    
    % SAVE uVECTOR DATA TO .VTK FORMAT %
    if find(strcmp({Output_Input{:,1}},'save_uVec')) >= 0
        Output_Params(11) = Output_Input{find(strcmp({Output_Input{:,1}},'save_uVec')),2};
    else
        Output_Params(11) = 1;
    end
    
    % SAVE uMAG TO .VTK FORMAT %
    if find(strcmp({Output_Input{:,1}},'save_uMag')) >= 0
        Output_Params(12) = Output_Input{find(strcmp({Output_Input{:,1}},'save_uMag')),2};
    else
        Output_Params(12) = 1;
    end

    % SAVE uX SCALAR DATA TO .VTK FORMAT %
    if find(strcmp({Output_Input{:,1}},'save_uX')) >= 0
        Output_Params(13) = Output_Input{find(strcmp({Output_Input{:,1}},'save_uX')),2};
    else
        Output_Params(13) = 1;
    end
    
    % SAVE uY SCALAR TO .VTK FORMAT %
    if find(strcmp({Output_Input{:,1}},'save_uY')) >= 0
        Output_Params(14) = Output_Input{find(strcmp({Output_Input{:,1}},'save_uY')),2};
    else
        Output_Params(14) = 1;
    end
    
    
    % SAVE fMAG TO .VTK FORMAT %
    if find(strcmp({Output_Input{:,1}},'save_fMag')) >= 0
        Output_Params(15) = Output_Input{find(strcmp({Output_Input{:,1}},'save_fMag')),2};
    else
        Output_Params(15) = 1;
    end

    % SAVE fX SCALAR DATA TO .VTK FORMAT %
    if find(strcmp({Output_Input{:,1}},'save_fX')) >= 0
        Output_Params(16) = Output_Input{find(strcmp({Output_Input{:,1}},'save_fX')),2};
    else
        Output_Params(16) = 1;
    end
    
    % SAVE fY SCALAR TO .VTK FORMAT %
    if find(strcmp({Output_Input{:,1}},'save_fY')) >= 0
        Output_Params(17) = Output_Input{find(strcmp({Output_Input{:,1}},'save_fY')),2};
    else
        Output_Params(17) = 1;
    end 
    
    % SAVE LAG STRUCTURE SCALAR DATA TO .VTK FORMAT %
    if find(strcmp({Output_Input{:,1}},'save_hier')) >= 0
        Output_Params(18) = Output_Input{find(strcmp({Output_Input{:,1}},'save_hier')),2};
    else
        Output_Params(18) = 1;
    end  
    
    
catch
    fprintf('\n \n');
    error('OUTPUT Information Improperly Declared in input2d file');  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Initializes zero Lagrangian structure fiber models and zero for
%           output information
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Lag_Struct_Params,Output_Params,Con_Params] = please_Initialize_Off()

Lag_Struct_Params(1) = 0;         % Springs: 0 (for no) or 1 (for yes) 
Lag_Struct_Params(2) = 0;         % Update_Springs: 0 (for no) or 1 (for yes)
Lag_Struct_Params(3) = 0;         % Target_Pts: 0 (for no) or 1 (for yes)
Lag_Struct_Params(4) = 0;         % Update_Target_Pts: 0 (for no) or 1 (for yes)
Lag_Struct_Params(5) = 0;         % Beams (Torsional Springs): 0 (for no) or 1 (for yes)
Lag_Struct_Params(6) = 0;         % Update_Beams (Torsional Springs): 0 (for no) or 1 (for yes)
Lag_Struct_Params(7) = 0;         % Non-Invariant Beams: 0 (for no) or 1 (for yes)
Lag_Struct_Params(8) = 0;         % Update_NonInvariant_Beams: 0 (for no) or 1 (for yes)
Lag_Struct_Params(9) = 0;         % Muscle Activation (Length/Tension-Hill Model): 0 (for no) or 1 (for yes)
Lag_Struct_Params(10) = 0;        % Muscle Activation 3-ELEMENT HILL MODEL w/ Length-Tension/Force-Velocity: 0 (for no) or 1 (for yes)
Lag_Struct_Params(11) = 0;        % Arbirtary External Force Onto Fluid Grid: 0 (for no) or 1 (for yes)
Lag_Struct_Params(12) = 0;        % Tracer Particles: 0 (for no) or 1 (for yes)
Lag_Struct_Params(13)= 0;         % Mass Points: 0 (for no) or 1 (for yes)
Lag_Struct_Params(14)= 0;         % Gravity: 0 (for no) or 1 (for yes)
Lag_Struct_Params(15)= 0;         % x-Component of Gravity vector
Lag_Struct_Params(16)= 0;         % y-Component of Gravity Vector
Lag_Struct_Params(17)= 0;         % Porous Media: 0 (for no) or 1 (for yes)
Lag_Struct_Params(18)= 0;         % Electrophysiology Model (FitzHugh-Nagumo)
Lag_Struct_Params(19)= 0;         % Damped Springs: 0 (for no) or 1 (for yes)
Lag_Struct_Params(20)= 0;         % Update_Damped_Springs: 0 (for no) or 1 (for yes)
Lag_Struct_Params(21)= 0;         % Boussinesq: 0 (for no) or 1 (for yes)
Lag_Struct_Params(22)= 0;         % expansion coefficient for Boussinesq approx.
Lag_Struct_Params(23)= 0;         % user-defined general force model: 0 (for no) or 1 (for yes)
Lag_Struct_Params(24)= 0;         % poroelastic media: 0 (for no) or 1 (for yes)
Lag_Struct_Params(25)= 0;         % coagulation model: 0 (for no) or 1 (for yes)

Con_Params(1)= 0;         % Background Concentration Gradient: 0 (for no) or 1 (for yes)
Con_Params(2)= 0;         % Diffusion Coefficient
Con_Params(3)= 0;         % Advection scheme: 0 (for 1st O Upwind) or 1 (for 3rd O WENO))
Con_Params(4)= 0;         % Concentration source: 0 (for no) or 1 (for constant) or 2 (for limited))
Con_Params(5)= 0;         % Concentration source rate
Con_Params(6)= 0;         % Concentration saturation limit
%Con_Params(7)= 0;         % Concentration periodic boundary: 0 (for no) or 1 (for yes) 


Output_Params(1) = 1000;          % Print Dump (How often to plot)
Output_Params(2) = 0;             % Plot in Matlab? (1=YES,0=NO) 
Output_Params(3) = 0;             % Plot LAGRANGIAN PTs ONLY in Matlab
Output_Params(4) = 0;             % Plot LAGRANGIAN PTs + VELOCITY FIELD in Matlab
Output_Params(5) = 0;             % Plot LAGRANGIAN PTs + VORTICITY colormap in Matlab
Output_Params(6) = 0;             % Plot LAGRANGIAN PTs + MAGNITUDE OF VELOCITY colormap in Matlab
Output_Params(7) = 0;             % Plot LAGRANGIAN PTs + PRESSURE colormap in Matlab
Output_Params(8) = 0;             % Plot LAGRANGIAN PTs + CONCENTRATION colormap in Matlab



