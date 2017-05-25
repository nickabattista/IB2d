#-------------------------------------------------------------------------------------------------------------------#
#
# IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
# 	fluid-structure interaction models. This version of the code is based off of
#	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
#
# Author: Nicholas A. Battista
# Email:  nick.battista@unc.edu
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

from please_Read_input2d_File import *
import numpy as np
import sys

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

    # GIVEN IN MATLAB INDEXING NOTATION:
    # Output_Params(1): print_dump
    #              (2): plot_Matlab
    #              (3): plot_LagPts
    #              (4): plot_Velocity
    #              (5): plot_Vorticity
    #              (6): plot_MagVelocity
    #              (7): plot_Pressure
    #              (8):  save_Vorticity 
    #              (9):  save_Pressure 
    #              (10): save_uVec 
    #              (11): save_uMag 
    #              (12): save_uX 
    #              (13): save_uY 
    #              (14): save_fMag 
    #              (15): save_fX 
    #              (16): save_fY 
    #              (17): save_hier 

    # Initialize 
    Output_Params = np.zeros(17)

    try: 

        # PRINT_DUMP
        try:
            ind = Output_Input[0][:].index('print_dump')
            Output_Params[0] = Output_Input[1][ind] # MATLAB: Output_Input{find(strcmp({Output_Input{:,1}},'print_dump ')),2};
        except ValueError:
            try:
                ind = Output_Input[0][:].index('print_dum')
                Output_Params[0] = Output_Input[1][ind] # MATLAB: Output_Input{find(strcmp({Output_Input{:,1}},'print_dump ')),2};
            except ValueError:
                Output_Params[0] = 100  # DEFAULT

        
        #%%%%%                                         %%%%%%
        #%%%%% OUTPUT INFO PLOTTING DIRECTLY IN MATLAB %%%%%%
        #%%%%%                                         %%%%%%


        # plot_MATLAB flag
        try:
            ind = Output_Input[0][:].index('plot_Matlab')
            Output_Params[1] = Output_Input[1][ind] # MATLAB: Output_Input{find(strcmp({Output_Input{:,1}},'plot_Matlab ')),2};
        except ValueError:
            try:
                ind = Output_Input[0][:].index('plot_Matla')
                Output_Params[1] = Output_Input[1][ind] # MATLAB: Output_Input{find(strcmp({Output_Input{:,1}},'plot_Matlab ')),2};
            except ValueError:
                Output_Params[1] = 0    # DEFAULT
    


        # plot_LagPts flag
        try:
            ind = Output_Input[0][:].index('plot_LagPts')
            Output_Params[2] = Output_Input[1][ind] # MATLAB: Output_Input{find(strcmp({Output_Input{:,1}},'plot_Matlab ')),2};
        except ValueError:
            try:
                ind = Output_Input[0][:].index('plot_LagPt')
                Output_Params[2] = Output_Input[1][ind] # MATLAB: Output_Input{find(strcmp({Output_Input{:,1}},'plot_Matlab ')),2};
            except ValueError:
                Output_Params[2] = 0    # DEFAULT


        # plot_Velocity flag
        try:
            ind = Output_Input[0][:].index('plot_Velocity')
            Output_Params[3] = Output_Input[1][ind] # MATLAB: Output_Input{find(strcmp({Output_Input{:,1}},'plot_Matlab ')),2};
        except ValueError:
            try:
                ind = Output_Input[0][:].index('plot_Velocit')
                Output_Params[3] = Output_Input[1][ind] # MATLAB: Output_Input{find(strcmp({Output_Input{:,1}},'plot_Matlab ')),2};
            except ValueError:
                Output_Params[3] = 0    # DEFAULT


        # plot_Vorticity flag
        try:
            ind = Output_Input[0][:].index('plot_Vorticity')
            Output_Params[4] = Output_Input[1][ind] # MATLAB: Output_Input{find(strcmp({Output_Input{:,1}},'plot_Matlab ')),2};
        except ValueError:
            try:
                ind = Output_Input[0][:].index('plot_Vorticit')
                Output_Params[4] = Output_Input[1][ind] # MATLAB: Output_Input{find(strcmp({Output_Input{:,1}},'plot_Matlab ')),2};
            except ValueError:
                Output_Params[4] = 0    # DEFAULT


        # plot_MagVelocity flag
        try:
            ind = Output_Input[0][:].index('plot_MagVelocity')
            Output_Params[5] = Output_Input[1][ind] # MATLAB: Output_Input{find(strcmp({Output_Input{:,1}},'plot_Matlab ')),2};
        except ValueError:
            try:
                ind = Output_Input[0][:].index('plot_MagVelocit')
                Output_Params[5] = Output_Input[1][ind] # MATLAB: Output_Input{find(strcmp({Output_Input{:,1}},'plot_Matlab ')),2};
            except ValueError:
                Output_Params[5] = 0    # DEFAULT


        # plot_Pressure flag
        try:
            ind = Output_Input[0][:].index('plot_Pressure')
            Output_Params[6] = Output_Input[1][ind] # MATLAB: Output_Input{find(strcmp({Output_Input{:,1}},'plot_Matlab ')),2};
        except ValueError:
            try:
                ind = Output_Input[0][:].index('plot_Pressur')
                Output_Params[6] = Output_Input[1][ind] # MATLAB: Output_Input{find(strcmp({Output_Input{:,1}},'plot_Matlab ')),2};
            except ValueError:
                Output_Params[6] = 0    # DEFAULT

        #%%%%%                                                %%%%%%
        #%%%%% OUTPUT INFO FOR WHAT GETS SAVED TO .VTK FORMAT %%%%%%
        #%%%%%                                                %%%%%%

        # save_Vorticity flag
        try:
            ind = Output_Input[0][:].index('save_Vorticity')
            Output_Params[7] = Output_Input[1][ind] 
        except ValueError:
            try:
                ind = Output_Input[0][:].index('save_Vorticit')
                Output_Params[7] = Output_Input[1][ind] 
            except ValueError:
                Output_Params[7] = 1    # DEFAULTS TO ALWAYS SAVE DATA        

        # save_Pressure flag
        try:
            ind = Output_Input[0][:].index('save_Pressure')
            Output_Params[8] = Output_Input[1][ind] 
        except ValueError:
            try:
                ind = Output_Input[0][:].index('save_Pressur')
                Output_Params[8] = Output_Input[1][ind] 
            except ValueError:
                Output_Params[8] = 1    # DEFAULTS TO ALWAYS SAVE DATA  


        # save_uVec flag
        try:
            ind = Output_Input[0][:].index('save_uVec')
            Output_Params[9] = Output_Input[1][ind] 
        except ValueError:
            try:
                ind = Output_Input[0][:].index('save_uVe')
                Output_Params[9] = Output_Input[1][ind] 
            except ValueError:
                Output_Params[9] = 1    # DEFAULTS TO ALWAYS SAVE DATA        

        # save_uMag flag
        try:
            ind = Output_Input[0][:].index('save_uMag')
            Output_Params[10] = Output_Input[1][ind] 
        except ValueError:
            try:
                ind = Output_Input[0][:].index('save_uMa')
                Output_Params[10] = Output_Input[1][ind] 
            except ValueError:
                Output_Params[10] = 1    # DEFAULTS TO ALWAYS SAVE DATA                 

        # save_uX flag
        try:
            ind = Output_Input[0][:].index('save_uX')
            Output_Params[11] = Output_Input[1][ind] 
        except ValueError:
            try:
                ind = Output_Input[0][:].index('save_uX')
                Output_Params[11] = Output_Input[1][ind] 
            except ValueError:
                Output_Params[11] = 1    # DEFAULTS TO ALWAYS SAVE DATA        

        # save_uY flag
        try:
            ind = Output_Input[0][:].index('save_uY')
            Output_Params[12] = Output_Input[1][ind] 
        except ValueError:
            try:
                ind = Output_Input[0][:].index('save_uY')
                Output_Params[12] = Output_Input[1][ind] 
            except ValueError:
                Output_Params[12] = 1    # DEFAULTS TO ALWAYS SAVE DATA 


        # save_fMag flag
        try:
            ind = Output_Input[0][:].index('save_fMag')
            Output_Params[13] = Output_Input[1][ind] 
        except ValueError:
            try:
                ind = Output_Input[0][:].index('save_fMa')
                Output_Params[13] = Output_Input[1][ind] 
            except ValueError:
                Output_Params[13] = 1    # DEFAULTS TO ALWAYS SAVE DATA                 

        # save_fX flag
        try:
            ind = Output_Input[0][:].index('save_fX')
            Output_Params[14] = Output_Input[1][ind] 
        except ValueError:
            try:
                ind = Output_Input[0][:].index('save_fX')
                Output_Params[14] = Output_Input[1][ind] 
            except ValueError:
                Output_Params[14] = 1    # DEFAULTS TO ALWAYS SAVE DATA   

         # save_fY flag
        try:
            ind = Output_Input[0][:].index('save_fY')
            Output_Params[15] = Output_Input[1][ind] 
        except ValueError:
            try:
                ind = Output_Input[0][:].index('save_fY')
                Output_Params[15] = Output_Input[1][ind] 
            except ValueError:
                Output_Params[15] = 1    # DEFAULTS TO ALWAYS SAVE DATA                

         # save_hier flag
        try:
            ind = Output_Input[0][:].index('save_hier')
            Output_Params[16] = Output_Input[1][ind] 
        except ValueError:
            try:
                ind = Output_Input[0][:].index('save_hier')
                Output_Params[16] = Output_Input[1][ind] 
            except ValueError:
                Output_Params[16] = 1    # DEFAULTS TO ALWAYS SAVE DATA




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

# GIVEN IN MATLAB INDEXING NOTATION:
# Lag_Struct_Params(1): springs
#                  (2): update_springs
#                  (3): target points
#                  (4): update_target_points
#                  (5): beams (torsional beams)
#                  (6): update_beams
#                  (7): nonInvariant_beams
#                  (8): update_nonInv_beams
#                  (9): FV_LT_muscle
#                  (10): 3_element_muscle
#                  (11): arb_ext_force
#                  (12): tracers
#                  (13): mass_pts
#                  (14): gravity
#                  (15): x_gravity_vec_comp
#                  (16): y_gravity_vec_comp
#                  (17): porous_media
#                  (18): concentration
#                  (19): electro_phys
#                  (20): damped_springs
#                  (21): update_damp_springs
#                  (22): boussinesq
#                  (23): expansion_coeff
#                  (24): user_force_model
#                   .         .
#                   .         .
#                   .         .

    # Initialize 
    Lag_Struct_Params = np.zeros(25)

    try: 

        # SPRINGS
        try:
            ind = Lag_Struct_Input[0][:].index('springs')
            Lag_Struct_Params[0] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'springs ')),2};
        except ValueError:
            try:
                ind = Lag_Struct_Input[0][:].index('spring')
                Lag_Struct_Params[0] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'springs ')),2};
            except ValueError:
                Lag_Struct_Params[0] = 0

        # UPDATE_SPRINGS
        try:        
            ind = Lag_Struct_Input[0][:].index('update_springs')
            Lag_Struct_Params[1] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'update_springs ')),2};
        except ValueError:
            try:
                ind = Lag_Struct_Input[0][:].index('update_spring')
                Lag_Struct_Params[1] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'update_springs ')),2};
            except ValueError:
                Lag_Struct_Params[1] = 0

       # TARGET PTS
        try:         
            ind = Lag_Struct_Input[0][:].index('target_pts')
            Lag_Struct_Params[2] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'target_pts ')),2};
        except ValueError:
            try:
                ind = Lag_Struct_Input[0][:].index('target_pt')
                Lag_Struct_Params[2] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'target_pts ')),2};
            except ValueError:
                Lag_Struct_Params[2] = 0


        # UPDATE_TARGET PTS
        try:         
            ind = Lag_Struct_Input[0][:].index('update_target')
            Lag_Struct_Params[3] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'update_target_pts ')),2};
        except ValueError:
            try:
                ind = Lag_Struct_Input[0][:].index('update_targe')
                Lag_Struct_Params[3] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'update_target_pts ')),2};
            except ValueError:
                Lag_Struct_Params[3] = 0
    

        # BEAM (TORSIONAL SPRING)
        try:
            ind = Lag_Struct_Input[0][:].index('beams')
            Lag_Struct_Params[4] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'beams ')),2};
        except ValueError:
            try:
                ind = Lag_Struct_Input[0][:].index('beam')
                Lag_Struct_Params[4] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'beams ')),2};
            except ValueError:
                Lag_Struct_Params[4] = 0

        # UPDATE_BEAMS (TORSIONAL SPRINGS)
        try:
            ind = Lag_Struct_Input[0][:].index('update_beams')
            Lag_Struct_Params[5] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'update_beams ')),2};
        except ValueError:
            try:
                ind = Lag_Struct_Input[0][:].index('update_beam')
                Lag_Struct_Params[5] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'update_beams ')),2};
            except ValueError:
                Lag_Struct_Params[5] = 0

        # NON-INVARIANT BEAMS
        try:
            ind = Lag_Struct_Input[0][:].index('nonInvariant_beams')
            Lag_Struct_Params[6] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'nonInvariant_beams ')),2};
        except ValueError:
            try:
                ind = Lag_Struct_Input[0][:].index('nonInvariant_beam')
                Lag_Struct_Params[6] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'nonInvariant_beams ')),2};
            except ValueError:
                Lag_Struct_Params[6] = 0


        # UPDATE_NON-INVARIANT_BEAMS
        try:
            ind = Lag_Struct_Input[0][:].index('update_nonInv_beams')
            Lag_Struct_Params[7] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'update_nonInv_beams ')),2};
        except ValueError:
            try:
                ind = Lag_Struct_Input[0][:].index('update_nonInv_beam')
                Lag_Struct_Params[7] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'update_nonInv_beams ')),2};
            except ValueError:
                Lag_Struct_Params[7] = 0
    

        # FV-LT MUSCLE
        try:
            ind = Lag_Struct_Input[0][:].index('FV_LT_muscle')
            Lag_Struct_Params[8] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'FV_LT_muscle ')),2};
        except ValueError:
            try:
                ind = Lag_Struct_Input[0][:].index('FV_LT_muscl')
                Lag_Struct_Params[8] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'FV_LT_muscle ')),2};
            except ValueError:
                Lag_Struct_Params[8] = 0


        # 3-ELEMENT HILL MUSCLE
        try:
            ind = Lag_Struct_Input[0][:].index('3_element_muscle')
            Lag_Struct_Params[9] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'3_element_muscle ')),2};
        except ValueError:
            try:
                ind = Lag_Struct_Input[0][:].index('3_element_muscl')
                Lag_Struct_Params[9] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'3_element_muscle ')),2};
            except ValueError:
                Lag_Struct_Params[9] = 0


        # ARBITRARY_EXTERNAL_FORCE
        try:
            ind = Lag_Struct_Input[0][:].index('arb_ext_force')
            Lag_Struct_Params[10] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'arb_ext_force ')),2};
        except ValueError:
            try:
                ind = Lag_Struct_Input[0][:].index('arb_ext_forc')
                Lag_Struct_Params[10] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'arb_ext_force ')),2};
            except ValueError:
                Lag_Struct_Params[10] = 0
    

        # TRACERS
        try:
            ind = Lag_Struct_Input[0][:].index('tracers')
            Lag_Struct_Params[11] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'tracers ')),2};
        except ValueError:
            try:
                ind = Lag_Struct_Input[0][:].index('tracer')
                Lag_Struct_Params[11] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'tracers ')),2};
            except ValueError:
                Lag_Struct_Params[11] = 0

    
        # MASS_POINTS
        try:
            ind = Lag_Struct_Input[0][:].index('mass_pts')
            Lag_Struct_Params[12] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'mass_pts ')),2};
        except ValueError:
            try:
                ind = Lag_Struct_Input[0][:].index('mass_pt')
                Lag_Struct_Params[12] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'mass_pts ')),2};
            except ValueError:
                Lag_Struct_Params[12] = 0


        # GRAVITY
        try:
            ind = Lag_Struct_Input[0][:].index('gravity')
            Lag_Struct_Params[13] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'gravity ')),2};
        except ValueError:
            try:
                ind = Lag_Struct_Input[0][:].index('gravit')
                Lag_Struct_Params[13] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'gravity ')),2};
            except ValueError:
                Lag_Struct_Params[13] = 0


         # x-GRAVITY COMPONENT
        try:
            ind = Lag_Struct_Input[0][:].index('x_gravity_vec_comp')
            Lag_Struct_Params[14] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'x_gravity_vec_comp ')),2};
        except ValueError:
            try:
                ind = Lag_Struct_Input[0][:].index('x_gravity_vec_com')
                Lag_Struct_Params[14] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'x_gravity_vec_comp ')),2};
            except ValueError:
                Lag_Struct_Params[14] = 0   


        # y-GRAVITY COMPONENT
        try:
            ind = Lag_Struct_Input[0][:].index('y_gravity_vec_comp')
            Lag_Struct_Params[15] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'y_gravity_vec_comp ')),2};
        except ValueError:
            try:
                ind = Lag_Struct_Input[0][:].index('y_gravity_vec_com')
                Lag_Struct_Params[15] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'y_gravity_vec_comp ')),2};
            except ValueError:
                Lag_Struct_Params[15] = 0   

    
        # POROUS_MEDIA
        try:
            ind = Lag_Struct_Input[0][:].index('porous_media')
            Lag_Struct_Params[16] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'porous_media ')),2};
        except ValueError:
            try:
                ind = Lag_Struct_Input[0][:].index('porous_medi')
                Lag_Struct_Params[16] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'porous_media ')),2};
            except ValueError:
                Lag_Struct_Params[16] = 0  


        # CONCENTRATION
        try:
            ind = Lag_Struct_Input[0][:].index('concentration')
            Lag_Struct_Params[17] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'concentration ')),2};
        except ValueError:
            try:
                ind = Lag_Struct_Input[0][:].index('concentratio')
                Lag_Struct_Params[17] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'concentration ')),2};
            except ValueError:
                Lag_Struct_Params[17] = 0 
    

        # ELECTROPHYS
        try:
            ind = Lag_Struct_Input[0][:].index('electro_phys')
            Lag_Struct_Params[18] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'electro_phys ')),2};
        except ValueError:
            try:
                ind = Lag_Struct_Input[0][:].index('electro_phy')
                Lag_Struct_Params[18] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'electro_phys ')),2};
            except ValueError:
                Lag_Struct_Params[18] = 0 


        # DAMPED_SPRINGS
        try:
            ind = Lag_Struct_Input[0][:].index('damped_springs')
            Lag_Struct_Params[19] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'damped_springs ')),2};
        except ValueError:
            try:
                ind = Lag_Struct_Input[0][:].index('damped_spring')
                Lag_Struct_Params[19] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'damped_springs ')),2};
            except ValueError:
                Lag_Struct_Params[19] = 0 


        # UPDATE_DAMPED_SPRINGS
        try:
            ind = Lag_Struct_Input[0][:].index('update_damp_springs')
            Lag_Struct_Params[20] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'update_damp_springs ')),2};
        except ValueError:
            try:
                ind = Lag_Struct_Input[0][:].index('update_damp_spring')
                Lag_Struct_Params[20] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'update_damp_springs ')),2};
            except ValueError:
                Lag_Struct_Params[20] = 0    


        # BOUSSINESQ
        try:
            ind = Lag_Struct_Input[0][:].index('boussinesq')
            Lag_Struct_Params[21] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'boussinesq ')),2};
        except ValueError:
            try:
                ind = Lag_Struct_Input[0][:].index('boussines')
                Lag_Struct_Params[21] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'boussinesq ')),2};
            except ValueError:
                Lag_Struct_Params[21] = 0  


        # EXPANSION_COEFFICIENT_FOR_BOUSSINESQ
        try:
            ind = Lag_Struct_Input[0][:].index('expansion_coeff')
            Lag_Struct_Params[22] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'expansion_coeff ')),2};
        except ValueError:
            try:
                ind = Lag_Struct_Input[0][:].index('expansion_coef')
                Lag_Struct_Params[22] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'expansion_coeff ')),2};
            except ValueError:
                Lag_Struct_Params[22] = 0


        # EXPANSION_COEFFICIENT_FOR_BOUSSINESQ
        try:
            ind = Lag_Struct_Input[0][:].index('user_force_model')
            Lag_Struct_Params[23] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'user_force_model ')),2};
        except ValueError:
            try:
                ind = Lag_Struct_Input[0][:].index('user_force_mode')
                Lag_Struct_Params[23] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'user_force_model ')),2};
            except ValueError:
                Lag_Struct_Params[23] = 0

        # POROELASTIC_MEDIA
        try:
            ind = Lag_Struct_Input[0][:].index('poroelastic')
            Lag_Struct_Params[24] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'user_force_model ')),2};
        except ValueError:
            try:
                ind = Lag_Struct_Input[0][:].index('poroelasti')
                Lag_Struct_Params[24] = Lag_Struct_Input[1][ind] # MATLAB: Lag_Struct_Input{find(strcmp({Lag_Struct_Input{:,1}},'user_force_model ')),2};
            except ValueError:
                Lag_Struct_Params[24] = 0        

    



      
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
    