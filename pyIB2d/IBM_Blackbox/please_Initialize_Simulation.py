#-------------------------------------------------------------------------------------------------------------------#
#
# IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
# 	fluid-structure interaction models. This version of the code is based off of
#	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
#
# Author: Nicholas A. Battista
# Email:  nickabattista@gmail.com
# Date Created: May 27th, 2015
# Institution: UNC-CH
#
# This code is capable of creating Lagrangian Structures using:
# 	1. Springs
# 	2. Beams (*torsional springs)
# 	3. Target Points
#	4. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-Tension)")
#   5. etc, etc.
#
# One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting lengths, etc
# 
# There are a number of built in Examples, mostly used for teaching purposes. 
# 
# If you would like us to add a specific muscle model, please let Nick (nickabattista@gmail.edu) know.
#
#--------------------------------------------------------------------------------------------------------------------#

from please_Read_input2d_File import please_Read_input2d_File
import numpy as np
import sys
from collections import OrderedDict

##################################################################################
#
# FUNCTION: initializes FLUID parameters for IBM_Driver file
#
#################################################################################

def please_Initialize_Fluid_Inputs(Fluid_Input):

    # GIVEN IN MATLAB INDEXING NOTATION:
    # Fluid_Params[1]: mu
    #             [2]: density

    # Initialize 
    Fluid_Params = np.zeros(2)

    try: 
        ind = Fluid_Input[0][:].index('mu')
        Fluid_Params[0] = Fluid_Input[1][ind] # MATLAB: Fluid_Input{find(strcmp({Fluid_Input{:,1}},'mu ')),2};
    
        ind = Fluid_Input[0][:].index('rho')
        Fluid_Params[1] = Fluid_Input[1][ind] # MATLAB: Fluid_Input{find(strcmp({Fluid_Input{:,1}},'rho ')),2};

    except ValueError:

        print('\n\n * FLUID Parameters Improperly Declared in input2d file * \n\n')
        print(' * Check input2d file format: \n')
        print('      Fluid_Parameters { \n')
        print('      mu =  <dynamic viscosity (N*s/m^2)>\n')
        print('      rho = <density (kg/m^3)>\n')
        print('      }\n')
        print('ERROR ERROR ERROR\n\n')

        sys.exit(0)
    
    return Fluid_Params

##################################################################################
#
# FUNCTION: initializes GRID parameters for IBM_Driver file
#
#################################################################################

def please_Initialize_Grid_Inputs(Grid_Input):

    # GIVEN IN MATLAB INDEXING NOTATION:
    # Grid_Params(1): Nx
    #            (2): Ny
    #            (3): Lx
    #            (4): Ly
    #            (5): Supp

    # Initialize 
    Grid_Params = np.zeros(5)

    try: 
        ind = Grid_Input[0][:].index('Nx')
        Grid_Params[0] = Grid_Input[1][ind] # MATLAB: Grid_Input{find(strcmp({Grid_Input{:,1}},'Nx ')),2};
    
        ind = Grid_Input[0][:].index('Ny')
        Grid_Params[1] = Grid_Input[1][ind] # MATLAB: Grid_Input{find(strcmp({Grid_Input{:,1}},'Ny ')),2};

        ind = Grid_Input[0][:].index('Lx')
        Grid_Params[2] = Grid_Input[1][ind] # MATLAB: Grid_Input{find(strcmp({Grid_Input{:,1}},'Lx ')),2};
    
        ind = Grid_Input[0][:].index('Ly')
        Grid_Params[3] = Grid_Input[1][ind] # MATLAB: Grid_Input{find(strcmp({Grid_Input{:,1}},'Ly ')),2};

        ind = Grid_Input[0][:].index('supp')
        Grid_Params[4] = Grid_Input[1][ind] # MATLAB: Grid_Input{find(strcmp({Grid_Input{:,1}},'supp ')),2};

    except ValueError:      

        print('\n\n * GRID Parameters Improperly Declared in input2d file * \n\n')
        print(' * Check input2d file format: \n')
        print('      Grid_Parameters { \n')
        print('      Nx = <# of Eulerian Grid Pts. in x-Direction (MUST BE EVEN!)>\n')
        print('      Ny = <# of Eulerian Grid Pts. in y-Direction (MUST BE EVEN!)>\n')
        print('      Lx = <Length of Eulerian Grid in x-Direction (m)>\n')
        print('      Ly = <Length of Eulerian Grid in y-Direction (m)>\n')
        print('      supp =  <Choose dirac-delta support (KEEP IT EVEN! *only set up for supp = 4*)>\n')
        print('      }\n')
        print('ERROR ERROR ERROR\n\n')
        sys.exit(0)
    
    return Grid_Params   


##################################################################################
#
# FUNCTION: initializes TEMPORAL information for IBM_Driver file
#
#################################################################################

def please_Initialize_Time_Inputs(Time_Input):

    # GIVEN IN MATLAB INDEXING NOTATION:
    # Time_Params(1): Tfinal (end time of simulation)
    #            (2): dt (time-step)

    # Initialize 
    Time_Params = np.zeros(2)

    try: 
        ind = Time_Input[0][:].index('Tfinal')
        Time_Params[0] = Time_Input[1][ind] # MATLAB: Time_Input{find(strcmp({Time_Input{:,1}},'Tfinal ')),2};
    
        ind = Time_Input[0][:].index('dt')
        Time_Params[1] = Time_Input[1][ind] # MATLAB: Time_Input{find(strcmp({Time_Input{:,1}},'dt ')),2};

    except ValueError:

        print('\n\n * TEMPORAL Parameters Improperly Declared in input2d file * \n\n')
        print(' * Check input2d file format: \n')
        print('      Temporal_Information { \n')
        print('      Tfinal =  <final time of simulation>\n')
        print('      dt = <time-step value>\n')
        print('      }\n')
        print('ERROR ERROR ERROR\n\n')

        sys.exit(0)
    
    return Time_Params


##################################################################################
#
# FUNCTION: initializes OUTPUT information for IBM_Driver file
#
#################################################################################

def please_Initialize_Output_Inputs(Output_Input):

    #                                                                      #
    # SPECIFY HERE A DICTIONARY OF ALL PARAMETERS AND THEIR DEFAULT VALUE! #
    #                                                                      #
    PARAMS_DEFAULT = (('print_dump', 100), 
                      ('plot_Matlab', 0),
                      ('plot_LagPts', 0),
                      ('plot_Velocity', 0),
                      ('plot_Vorticity', 0),
                      ('plot_MagVelocity', 0),
                      ('plot_Pressure', 0),
                      ('save_Vorticity', 1),
                      ('save_Pressure', 1),
                      ('save_uVec', 1),
                      ('save_uMag', 1),
                      ('save_uX', 1),
                      ('save_uY', 1),
                      ('save_fMag', 1),
                      ('save_fX', 1),
                      ('save_fY', 1),
                      ('save_hier', 1)
    )

    PARAMS_DEFAULT = OrderedDict(PARAMS_DEFAULT)

    # Initialize 
    Output_Params = np.zeros(len(PARAMS_DEFAULT))

    try: 
        n = 0
        for param, default_val in PARAMS_DEFAULT.items():
            try:
                ind = Output_Input[0][:].index(param)
                Output_Params[n] = Output_Input[1][ind]
            except ValueError:
                try:
                    ind = Output_Input[0][:].index(param[:-1])
                    Output_Params[n] = Output_Input[1][ind]
                except ValueError:
                    Output_Params[n] = default_val
            n += 1

    except ValueError:
        
        print('\n\nERROR ERROR ERROR')
        print('\n\n * OUTPUT Parameters Improperly Declared in input2d file * \n\n')
        print(' * Check input2d file format: \n')
        print('      Output_Info { \n')
        print('      print_dump =  <# of time-steps btwn saving data>\n')
        print('      plot_Matlab = <0 or 1>\n')
        print('      plot_LagPts = <0 or 1>\n')
        print('      plot_Velocity = <0 or 1>\n')
        print('      plot_Vorticity = <0 or 1>\n')
        print('      plot_MagVelocity = <0 or 1>\n')
        print('      plot_Pressure = <0 or 1>\n')
        print('      }\n')
        print('ERROR ERROR ERROR\n\n')
        sys.exit(0)
    
    return Output_Params


##################################################################################
#
# FUNCTION: initializes Lag_Structure information for IBM_Driver file
#
#################################################################################

def please_Initialize_Lag_Structure_Inputs(Lag_Struct_Input):

    #                                                                      #
    # SPECIFY HERE A DICTIONARY OF ALL PARAMETERS AND THEIR DEFAULT VALUE! #
    #                                                                      #

    PARAMS_DEFAULT = (('springs', 0),
                      ('update_springs', 0),
                      ('target_pts', 0),
                      ('update_target', 0),
                      ('beams', 0),
                      ('update_beams', 0),
                      ('nonInvariant_beams', 0),
                      ('update_nonInv_beams', 0),
                      ('FV_LT_muscle', 0),
                      ('3_element_muscle', 0),
                      ('arb_ext_force', 0),
                      ('tracers', 0),
                      ('mass_pts', 0),
                      ('gravity', 0),
                      ('x_gravity_vec_comp', 0),
                      ('y_gravity_vec_comp', 0),
                      ('porous_media', 0),
                      ('concentration', 0),
                      ('electro_phys', 0),
                      ('damped_springs', 0),
                      ('update_damp_springs', 0),
                      ('boussinesq', 0),
                      ('expansion_coeff', 0),
                      ('user_force_model', 0),
                      ('poroelastic', 0),
                      ('brinkman', 0)
                      )

    PARAMS_DEFAULT = OrderedDict(PARAMS_DEFAULT)

    # Initialize 
    Lag_Struct_Params = np.zeros(len(PARAMS_DEFAULT))

    try: 
        n = 0
        for param, default_val in PARAMS_DEFAULT.items():
            try:
                ind = Lag_Struct_Input[0][:].index(param)
                Lag_Struct_Params[n] = Lag_Struct_Input[1][ind]
            except ValueError:
                try:
                    ind = Lag_Struct_Input[0][:].index(param[:-1])
                    Lag_Struct_Params[n] = Lag_Struct_Input[1][ind]
                except ValueError:
                    Lag_Struct_Params[n] = default_val
            n += 1
      
    except ValueError:
        
        print('\n\nERROR ERROR ERROR')
        print('\n\n * LAGRANGIAN STRUCTURE Selections Improperly Declared in input2d file * \n\n')
        print(' * Check input2d file format: \n')
        print('      Lag_Structure_Info { \n')
        print('      springs =  <0 or 1>\n')
        print('      update_springs = <0 or 1>\n')
        print('      target_pts = <0 or 1>\n')
        print('      update_target = <0 or 1>\n')
        print('      beams = <0 or 1>\n')
        print('      update_beams = <0 or 1>\n')
        print('      nonInvariant_beams = <0 or 1>\n')
        print('      update_nonInv_beams = <0 or 1>\n')
        print('      FV_LT_muscle =  <0 or 1>\n')
        print('      3_element_muscle = <0 or 1>\n')
        print('      arb_ext_force = <0 or 1>\n')
        print('      tracers = <0 or 1>\n')
        print('      mass_pts = <0 or 1>\n')
        print('      gravity = <0 or 1>\n')
        print('      x_gravity_vec_comp = <# (x-GRAVITY VECTOR COMPONENT) >\n')
        print('      y_gravity_vec_comp = <# (y-GRAVITY VECTOR COMPONENT) >\n')
        print('      porous_media =  <0 or 1>\n')
        print('      concentration = <0 or 1>\n')
        print('      electro_phys = <0 or 1>\n')
        print('      damped_springs = <0 or 1>\n')
        print('      update_damp_springs = <0 or 1>\n')
        print('      boussinesq = <0 or 1>\n')
        print('      expansion_coeff = <# (EXPANSION COEFFICIENT FOR BOUSSINESQ) >\n')
        print('      user_force_model = <0 or 1>\n')
        print('      poroelastic = <0 or 1>\n')
        print('      brinkman = <0 or 1>\n')
        print('      }\n')
        print('ERROR ERROR ERROR\n\n')
        sys.exit(0)
    
    return Lag_Struct_Params



##################################################################################
#
# FUNCTION: initializes LAGRANGIAN STRUCTURE NAME information for IBM_Driver file
#
#################################################################################

def please_Initialize_Lag_Name_Inputs(Lag_Name_Input):

    # Lag_Name_Params: string_name

    try: 
        ind = Lag_Name_Input[0][:].index('string_name')
        Lag_Name_Params = Lag_Name_Input[1][ind] # MATLAB: Lag_Name_Input{find(strcmp({Lag_Name_Input{:,1}},'string_name ')),2}
    
    except ValueError:

        print('\n\n * LAGRANGIAN NAME Parameter Improperly Declared in input2d file * \n\n')
        print(' * Check input2d file format: \n')
        print('      Lag_Name { \n')
        print('      string_name = <"Lagrangian_Structure_Name" (in quotation marks)>\n')
        print('      }\n')
        print('ERROR ERROR ERROR\n\n')

        sys.exit(0)
    
    return Lag_Name_Params


##################################################################################
#
# FUNCTION: reads in input2d files and initializes the simulation.
#
#################################################################################

def please_Initialize_Simulation():

    
    #
    # READ IN ALL INPUTS INTO CELLS FROM INPUT2D #
    #
    params = please_Read_input2d_File('input2d')

    #
    # EXTRACT INDIVIDUAL CELL GROUPS #
    #
    ind = params[0][:].index('Fluid_Parameters')
    Fluid_Input = params[1][ind][:]             # MATLAB: Fluid_Input = params{find(strcmp({params{:,1}},'Fluid_Parameters')),2}
    
    ind = params[0][:].index('Grid_Parameters')
    Grid_Input = params[1][ind][:]              # MATLAB: Grid_Input = params{find(strcmp({params{:,1}},'Grid_Parameters')),2}

    ind = params[0][:].index('Temporal_Information')
    Time_Input = params[1][ind][:]              # MATLAB: Time_Input = params{find(strcmp({params{:,1}},'Temporal_Information')),2}
    
    ind = params[0][:].index('Lag_Structure_Info')
    Lag_Struct_Input = params[1][ind][:]        # MATLAB: Lag_Struct_Input = params{find(strcmp({params{:,1}},'Lag_Structure_Info')),2}
    
    ind = params[0][:].index('Output_Info')
    Output_Input = params[1][ind][:]             # MATLAB: Output_Input = params{find(strcmp({params{:,1}},'Output_Info')),2}
    
    ind = params[0][:].index('Lag_Name')
    Lag_Name_Input = params[1][ind][:]          #MATLAB: Lag_Name_Input = params{find(strcmp({params{:,1}},'Lag_Name')),2}


    #
    # TESTING READ INPUTS
    #
    #print(Fluid_Input)
    #print(Grid_Input)
    #print(Time_Input)
    #print(Lag_Struct_Input)
    #print(Output_Info)
    #print(Lag_Name_Input)


    #
    # INITIALIZE PARAMETERS FOR IBM_DRIVER FILE #
    #
    Fluid_Params = please_Initialize_Fluid_Inputs(Fluid_Input)
    Grid_Params = please_Initialize_Grid_Inputs(Grid_Input)
    Time_Params = please_Initialize_Time_Inputs(Time_Input)
    Lag_Struct_Params = please_Initialize_Lag_Structure_Inputs(Lag_Struct_Input)
    Output_Params = please_Initialize_Output_Inputs(Output_Input)
    Lag_Name_Params = please_Initialize_Lag_Name_Inputs(Lag_Name_Input)
    
    #
    # TESTING READ INPUTS STORAGE VALUES
    #
    #print(Fluid_Params)
    #print(Grid_Params)
    #print(Time_Params)
    #print(Output_Params)
    #print(Lag_Struct_Params)
    #print(Lag_Name_Params)

    return Fluid_Params, Grid_Params, Time_Params, Lag_Struct_Params, Output_Params, Lag_Name_Params

if __name__ == "__main__":
    please_Initialize_Simulation()
    