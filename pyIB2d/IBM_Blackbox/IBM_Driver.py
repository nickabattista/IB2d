'''-------------------------------------------------------------------------

 IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
 	fluid-structure interaction models. This version of the code is based off of
	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.

 Author: Nicholas A. Battista
 Email:  nickabattista@gmail.com
 Date Created: May 27th, 2015
 Initial Python 3.5 port and VTK writes by Christopher Strickland
 Institution: UNC-CH

 This code is capable of creating Lagrangian Structures using:
 	1. Springs
 	2. Beams (*torsional springs)
 	3. Target Points
    4. Mass Points
    5. Porous Points
	6. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-Tension)")
    7. 3-Element Hill Muscle Model

 One is able to update those Lagrangian Structure parameters, e.g., 
 spring constants, resting lengths, etc
 
 There are a number of built in Examples, mostly used for teaching purposes. 
 
 If you would like us to add a specific muscle model, 
 please let Nick (nickabattista@gmail.com) know.

----------------------------------------------------------------------------'''
import numpy as np
from math import sqrt
import os
from Supp import *
from please_Find_Lagrangian_Forces_On_Eulerian_grid import\
    please_Find_Lagrangian_Forces_On_Eulerian_grid
from please_Update_Fluid_Velocity import please_Update_Fluid_Velocity
from please_Compute_Porous_Slip_Velocity import\
    please_Compute_Porous_Slip_Velocity
from please_Plot_Results import please_Plot_Results
from please_Compute_Normal_Tangential_Forces_On_Lag_Pts import \
    please_Compute_Normal_Tangential_Forces_On_Lag_Pts


# IF BOUSSINESQ
try:
    from please_Form_Boussinesq_Forcing_Terms import *
except ImportError:
    pass

#Here is the try import C part
try:
    import write
    C_flag = True
    print('Running with compiled C I/O library.')
except:
    C_flag = False

#Switch for vtk library writes (this will be overridden by C_flag)
try:
    import vtk
    from vtk.util import numpy_support
    vtk_lib_flag = True
except:
    vtk_lib_flag = False
    if not C_flag:
        print('Running without optimized IO libraries (VTK or C).')

###############################################################################
#
# FUNCTION: Actual DRIVER of the code, where the time-stepping occurs ->
#           gets called by main2d to do the "magic" :)
#
###############################################################################

def main(Fluid_Params,Grid_Params,Time_Params,Lag_Struct_Params,Output_Params,Lag_Name_Params):

    ''' 2D IMMERSED BOUNDARY SOLVER ON RECTANGULAR DOMAIN w/ PERIODIC BOUNDARIES

    Args:
        struct_name: structure name
        mu: dynamic viscosity
        rho: density
        grid_Info: dict of grid properties
        dt: time-step
        T_FINAL: final simulation time
        model_Info: dict of model structure properties
        
    Returns:
        X: here
        Y: are
        U: x-Eulerian grid velocity
        V: y-Eulerian grid velocity
        xLag: what do
        yLag: they do?
        

    x-Momentum Conservation: rho*u_t = -rho*u*u_x + rho*v*u_y + mu*laplacian(u) 
                                       - p_x + F_x
    y-Momentum Convervation: rho*v_t = -rho*u*v_x + rho*v*v_y + mu*laplacian(v) 
                                       - p_y + F_y
                                       
    Incompressibility: u_x + v_y = 0

    LagPts/dt = int{ u(x,t) delta( x - LagPts(s,t) ) dx }
    F_x = int{ fx(s,t) delta(x - LagPts(s,t)) ds }
    F_y = int{ fy(s,t) delta(x - LagPts(s,t)) ds }'''

    print('\n________________________________________________________________________________\n\n')
    print('\n---------------->>                 IB2d                      <<----------------\n')
    print('\n________________________________________________________________________________\n\n')
    print('If using the code for research purposes please cite the following two papers: \n')
    print('     [1] N.A. Battista, A.J. Baird, L.A. Miller, A mathematical model and MATLAB code for muscle-fluid-structure simulations, Integ. Comp. Biol. 55(5):901-11 (2015)\n')
    print('     [2] N.A. Battista, W.C. Strickland, L.A. Miller, IB2d a Python and MATLAB implementation of the immersed boundary method, Bioinspir. Biomim. 12(3):036003 (2017)\n')
    print('     [3] N.A. Battista, W.C. Strickland, A. Barrett, L.A. Miller, IB2d Reloaded: A more powerful Python and MATLAB implementation of the immersed boundary method. Math Meth Appl Sci. 1?26 (2017).')
    print('\n________________________________________________________________________________')
    print('\n\nNOTE: If running pyIB2d with Anaconda, please note that it will run faster with vtk installed.\n')
    print('          To install, type the following line into your terminal: <conda install -c menpo vtk> in your terminal')
    print('\n________________________________________________________________________________')
    print('\n\n\n |****** Prepping Immersed Boundary Simulation ******|\n')
    print('\n\n--> Reading input data for simulation...\n\n')
    
    #
    # ** IBM_DRIVER INPUT DEFINITIONS ** :
    #
    #               Fluid_Params[1]: mu
    #                           [2]: density
    #
    #               Grid_Params[1]: Nx
    #                          [2]: Ny
    #                          [3]: Lx
    #                          [4]: Ly
    #                          [5]: Supp
    #
    #               Time_Params[1]: Tfinal [end time of simulation]
    #                          [2]: dt [time-step]
    #
    #               Lag_Struct_Params[1]: springs
    #                                [2]: update_springs
    #                                [3]: target points
    #                                [4]: update_target_points
    #                                [5]: beams [torsional beams]
    #                                [6]: update_beams
    #                                 .         .
    #                                 .         .
    #                                 .         .
    #
    #               Output_Params[0]: print_dump
    #                            [1]: plot_Matlab
    #                            [2]: plot_LagPts
    #                            [3]: plot_Velocity
    #                            [4]: plot_Vorticity
    #                            [5]: plot_MagVelocity
    #                            [6]: plot_Pressure
    #                            [7]:  save_Vorticity 
    #                            [8]:  save_Pressure 
    #                            [9]: save_uVec 
    #                            [10]: save_uMag 
    #                            [11]: save_uX 
    #                            [12]: save_uY 
    #                            [13]: save_fMag 
    #                            [14]: save_fX 
    #                            [15]: save_fY 
    #                            [16]: save_hier 


    # SIMULATION NAME STRING TO RUN .vertex, .spring, etc. #
    struct_name = Lag_Name_Params

    # FLUID PARAMETER VALUES STORED #
    mu = Fluid_Params[0]      # Dynamic Viscosity
    rho = Fluid_Params[1]     # Density

    # TEMPORAL INFORMATION VALUES STORED #
    T_FINAL = Time_Params[0]     # Final simulation time
    dt = Time_Params[1]          # Time-step
    #NTime = floor[T_FINAL/dt]+1 # # of total time-steps [floor'd so exact number of time-steps]
    #dt = T_FINAL/NTime          # revised time-step [slightly perturbed dt, so exact # of time-steps are used]
    current_time = 0.0           # initialize start of simulation to time, 0 


    # GRID INFO #
    Nx =   int(Grid_Params[0])               # # of Eulerian pts. in x-direction
    Ny =   int(Grid_Params[1])                # # of Eulerian pts. in y-direction
    Lx =   Grid_Params[2]                # Length of Eulerian grid in x-coordinate
    Ly =   Grid_Params[3]                # Length of Eulerian grid in y-coordinate
    dx =   Grid_Params[2]/Grid_Params[0] # Spatial-size in x
    dy =   Grid_Params[3]/Grid_Params[1] # Spatial-size in y
    supp = Grid_Params[4]                # Delta-function support
    grid_Info = [Nx,Ny,Lx,Ly,dx,dy,supp] # Store for passing values into IB time-loop sub-functions
                                         # NOTE: grid_Info[8] = Nb, grid_Info[9] = ds [THEY GET STORED LATER]

    # PRINTING/PLOTTING INFO #
    pDump = Output_Params[0]          # Print [Plot] Dump interval
    pMatplotlib = Output_Params[1]    # Plot in Matlab? [1=YES,0=NO]
    lagPlot = Output_Params[2]        # Plot LAGRANGIAN PTs ONLY in Matlab
    velPlot = Output_Params[3]        # Plot LAGRANGIAN PTs + VELOCITY FIELD in Matlab
    vortPlot = Output_Params[4]       # Plot LAGRANGIAN PTs + VORTICITY colormap in Matlab
    uMagPlot = Output_Params[5]       # Plot LAGRANGIAN PTs + MAGNITUDE OF VELOCITY colormap in Matlab
    pressPlot = Output_Params[6]      # Plot LAGRANGIAN PTs + PRESSURE colormap in Matlab


    # MODEL STRUCTURE DATA STORED #
    springs_Yes = Lag_Struct_Params[0]           # Springs: 0 [for no] or 1 [for yes] 
    update_Springs_Flag = Lag_Struct_Params[1]   # Update_Springs: 0 [for no] or 1 [for yes]
    target_pts_Yes = Lag_Struct_Params[2]        # Target_Pts: 0 [for no] or 1 [for yes]
    update_Target_Pts = Lag_Struct_Params[3]     # Update_Target_Pts: 0 [for no] or 1 [for yes]
    beams_Yes = Lag_Struct_Params[4]             # Beams: 0 [for no] or 1 [for yes]
    update_Beams_Flag = Lag_Struct_Params[5]     # Update_Beams: 0 [for no] or 1 [for yes]
    nonInv_beams_Yes = Lag_Struct_Params[6]      # Beams [non-invariant]: 0 [for no] or 1 [for yes]
    update_nonInv_Beams_Flag = Lag_Struct_Params[7] # Update_nonInv_Beams: 0 [for no] or 1 [for yes]
    muscles_Yes = Lag_Struct_Params[8]           # FV-LT Muscles: 0 [for no] or 1 [for yes]
    hill_3_muscles_Yes = Lag_Struct_Params[9]    # Hill 3-Element Muscle: 0 [for no] or 1 [for yes]
    arb_ext_force_Yes = Lag_Struct_Params[10]     # Arbitrary External Force: 0 [for no] or 1 [for yes]
    tracers_Yes = Lag_Struct_Params[11]          # Tracers: 0 [for no] or 1 [for yes]
    mass_Yes = Lag_Struct_Params[12]             # Mass Points: 0 [for no] or 1 [for yes]
    gravity_Yes = Lag_Struct_Params[13]          # Gravity: 0 [for no] or 1 [for yes]
    # NOTE: Lag_Struct_Params[14],[15]:            <- components of gravity vector [if gravity, initialize them below]
    porous_Yes = Lag_Struct_Params[16]           # Porous Media: 0 [for no] or 1 [for yes]
    concentration_Yes = Lag_Struct_Params[17]    # Background Concentration Gradient: 0 [for no] or 1 [for yes]
    electro_phys_Yes = Lag_Struct_Params[18]     # Electrophysiology [FitzHugh-Nagumo]: 0 [for no] or 1 [for yes]
    d_Springs_Yes = Lag_Struct_Params[19]        # Damped Springs: 0 [for no] or 1 [for yes]
    update_D_Springs_Flag = Lag_Struct_Params[20]# Update_Damped_Springs: # 0 [for no] or 1 [for yes]
    boussinesq_Yes = Lag_Struct_Params[21]       # Boussinesq Approx.: 0 [for no] or 1 [for yes]
    exp_Coeff = Lag_Struct_Params[22]            # Expansion Coefficient [e.g., thermal, etc] for Boussinesq approx.
    general_force_Yes = Lag_Struct_Params[23]    # General User-Defined Force Term: 0 [for no] or 1 [for yes]
    poroelastic_Yes = Lag_Struct_Params[24]      # Poroelastic media 0 [for no] or 1 [for yes]
    brinkman_Yes = Lag_Struct_Params[25]         # Brinkman term in Momentum equation 0 [for no] or 1 [for yes]

    print(Lag_Struct_Params)

    # CLEAR INPUT DATA #


    # Temporal Information
    #NTime = np.floor(T_FINAL/dt)+1 # number of total time-steps,
                                # (floored, so exact number of time-steps)
    #dt = T_FINAL/NTime #time-step (slightly perturbed dt, so exact number of 
                       #time-steps are used
    #current_time = 0.0
    


    

    
    
    # Create EULERIAN Mesh (these assume periodicity in x and y)
    x = np.arange(0,Lx,dx)
    y = np.arange(0,Ly,dy)
    # Create x-Mesh
    #X = np.empty((Nx,x.size))
    #for ii in range(Nx):
    #    X[ii,] = x
    # Create y-Mesh
    #Y = np.empty((y.size,Ny))
    #for ii in range(Ny):
    #    Y[:,ii] = y
    # MATLAB SYNTAX: [X,Y] = meshgrid(0:dx:Lx-dx,0:dy:Ly-dy)
    X,Y = np.meshgrid(x,y)
    # MATLAB SYNTAX: [idX,idY] = meshgrid(0:Nx-1,0:Ny-1)  <--- INITIALIZE FOR FLUID SOLVER FFT FUNCTION    
    idX,idY = np.meshgrid(np.arange(0,Nx,1),np.arange(0,Ny,1)) # <- INITIALIZES FOR FLUID SOLVER FFT OPERATORS

    # # # # # HOPEFULLY WHERE I CAN READ IN INFO!!! # # # # #


    # READ IN LAGRANGIAN POINTS #
    Nb,xLag,yLag = read_Vertex_Points(struct_name)
    grid_Info.append(Nb)       # Total Number of Lagrangian Pts. (grid_Info[7]=Nb)
    xLag_P = xLag              # Initialize previous Lagrangian x-Values 
                               #   (for use in muscle-model)
    yLag_P = yLag              # Initialize previous Lagrangian y-Values 
                               #   (for use in muscle-model)

    #Lagrangian Structure Data
    ds = Lx/(2.*Nx)             # Lagrangian Spacing
    grid_Info.append(ds)        # grid_Info[8] = ds                  
    
    print('\n--> FIBER MODEL INCLUDES: \n')
                   
                            

       
       

        
        
    # READ IN SPRINGS (IF THERE ARE SPRINGS) #
    if springs_Yes:
        print('  - Springs and ...')
        if update_Springs_Flag == 0:
            print('                   NOT dynamically updating spring properties\n')
        else:
            print('                   dynamically updating spring properties\n')

        springs_info = read_Spring_Points(struct_name)
            #springs_info: col 1: starting spring pt (by lag. discretization)
            #              col 2: ending spring pt. (by lag. discretization)
            #              col 3: spring stiffness
            #              col 4: spring resting lengths
            #              col 5: degree of non-linearity (1=linear)
    else:
        springs_info = np.zeros((1,1))  #just to pass placeholder into 
            # "please_Find_Lagrangian_Forces_On_Eulerian_grid function"
    


    # READ IN BEAMS (IF THERE ARE TORSIONAL SPRINGS aka BEAMS) #
    if beams_Yes:
        print('  - Beams ("Torsional Springs") and ... ')
        if update_Beams_Flag == 0:
            print('                    NOT dynamically updating beam properties\n')
        else:
            print('                    dynamically updating beam properties\n')

        beams_info = read_Beam_Points(struct_name)
        #beams:      col 1: 1ST PT.
        #            col 2: MIDDLE PT. (where force is exerted)
        #            col 3: 3RD PT.
        #            col 4: beam stiffness
        #            col 5: curavture
    else:
        beams_info = 0



    # READ IN NON-INVARIANT BEAMS (IF THERE ARE NON-INVARIANT BEAMS) #
    if nonInv_beams_Yes:
        print('  - Beams ("Non-Invariant Beams") and ... ')
        if update_Beams_Flag == 0:
            print('                    NOT dynamically updating beam properties\n')
        else:
            print('                    dynamically updating beam properties\n')

        nonInv_beams_info = read_nonInv_Beam_Points(struct_name)
        #beams:      col 1: 1ST PT.
        #            col 2: MIDDLE PT. (where force is exerted)
        #            col 3: 3RD PT.
        #            col 4: beam stiffness
        #            col 5: x-curavture
        #            col 6: y-curvature
    else:
        nonInv_beams_info = 0    


    # READ IN DAMPED SPRINGS (IF THERE ARE DAMPED SPRINGS) #
    if d_Springs_Yes:
        print('  - Damped Springs and ...')
        if update_D_Springs_Flag == 0:
            print('                   NOT dynamically updating damped spring properties\n')
        else:
            print('                   dynamically updating spring properties\n')

        d_springs_info = read_Damped_Spring_Points(struct_name) 
            #springs_info: col 1: starting spring pt (by lag. discretization)
            #              col 2: ending spring pt. (by lag. discretization)
            #              col 3: spring stiffness
            #              col 4: spring resting lengths
            #              col 5: damping coefficient
    else:
        d_springs_info = np.zeros((1,1))  #just to pass placeholder into 
            # "please_Find_Lagrangian_Forces_On_Eulerian_grid function"
    
    

    # READ IN TARGET POINTS (IF THERE ARE TARGET PTS) #
    if target_pts_Yes:
        print('  - Target Pts. and ...')
        if update_Target_Pts == 0:
            print('                 NOT dynamically updating target point properties\n')
        else:
            print('                 dynamically updating target point properties\n')
        
        target_aux = read_Target_Points(struct_name)
        #target_aux: col 0: Lag Pt. ID w/ Associated Target Pt.
        #            col 1: target STIFFNESSES
        
        # initialize target_info
        target_info = np.empty((target_aux.shape[0],4))
        
        target_info[:,0] = target_aux[:,0] #Stores Lag-Pt IDs in col vector
        # Stores Original x-Lags and y-Lags as x/y-Target Pt. Identities
        target_info[:,1] = xLag[target_info[:,0].astype('int')]
        target_info[:,2] = yLag[target_info[:,0].astype('int')]
        
        target_info[:,3] = target_aux[:,1] #Stores Target Stiffnesses 
    else:
        target_info = np.zeros((1,1))
    
    

    
    
    
    # READ IN MASS POINTS (IF THERE ARE MASS PTS) #
    if mass_Yes:
        print('  - Mass Pts. with ')
        if gravity_Yes == 0:
            print('          NO artificial gravity\n')
        else:
            print('          artificial gravity\n')

        mass_aux = read_Mass_Points(struct_name)
        #target_aux: col 0: Lag Pt. ID w/ Associated Mass Pt.
        #            col 1: "Mass-spring" stiffness parameter
        #            col 2: "MASS" value parameter
        
        if ( mass_aux.size/3.0 < 2.0 ):  # checks if only 1 mass point

            # initialize mass_info
            mass_info = np.empty((1,5))
        
            mass_info[0,0] = mass_aux[0] #Stores Lag-Pt IDs in col vector
            
            #Stores Original x-Lags and y-Lags as x/y-Mass Pt. Identities
            if ( np.isscalar(xLag) ):
                mass_info[0,1] = xLag
                mass_info[0,2] = yLag
            else:
                mass_info[0,1] = xLag[mass_info[0].astype('int')]
                mass_info[0,2] = yLag[mass_info[0].astype('int')]

            mass_info[0,3] = mass_aux[1]   #Stores "mass-spring" parameter 
            mass_info[0,4] = mass_aux[2]   #Stores "MASS" value parameter
        
        else:
            # initialize mass_info
            mass_info = np.empty((mass_aux.shape[0],5))
        
            mass_info[:,0] = mass_aux[:,0] #Stores Lag-Pt IDs in col vector
            #Stores Original x-Lags and y-Lags as x/y-Mass Pt. Identities
            mass_info[:,1] = xLag[mass_info[:,0].astype('int')]
            mass_info[:,2] = yLag[mass_info[:,0].astype('int')]
        
            mass_info[:,3] = mass_aux[:,1]   #Stores "mass-spring" parameter 
            mass_info[:,4] = mass_aux[:,2]   #Stores "MASS" value parameter


    else:
        mass_info = np.zeros((1,1))






    
    # READ IN POROUS MEDIA INFO (IF THERE IS POROSITY) #
    if porous_Yes:
        print('  - Porous Points\n')
        porous_aux = read_Porous_Points(struct_name)
        #porous_aux: col 1: Lag Pt. ID w/ Associated Porous Pt.
        #            col 2: Porosity coefficient
        #            col 3: Porous Pt. 'Stencil' flag for derivative
        
        # initizlize porous_info
        porous_info = np.empty((porous_aux.shape[0],5))
        
        porous_info[:,0] = porous_aux[:,0] #Stores Lag-Pt IDs in col vector
        # Stores Original x-Lags and y-Lags as x/y-Porous Pt. Identities
        porous_info[:,1] = xLag[porous_info[:,0].astype('int')]
        porous_info[:,2] = yLag[porous_info[:,0].astype('int')]
        
        porous_info[:,3] = porous_aux[:,1] # Stores Porosity Coefficient 
        porous_info[:,4] = porous_aux[:,2] # Stores Porosity Stencil Flag
    else:
        porous_info = np.zeros((1,1))



    # READ IN PORO-ELASTIC MEDIA INFO (IF THERE IS PORO-ELASTICITY) #
    if ( poroelastic_Yes ):
        print('  -Poroelastic media\n')
        poroelastic_info = read_PoroElastic_Points(struct_name)
        F_Poro = np.zeros( ( len(poroelastic_info), 2) )   # Initialization
        #poroelastic_info: col 1: Lag Pt. ID w/ Associated Porous Pt.
        #                  col 2: Porosity coefficient
    else:
        poroelastic_info = np.zeros((1,1))
        F_Poro = np.zeros((1,1))

  


    # READ IN MUSCLES (IF THERE ARE MUSCLES) #
    if muscles_Yes:
        print('  - MUSCLE MODEL (Force-Velocity / Length-Tension Model)\n')
        muscles_info = read_Muscle_Points(struct_name)
            #         muscles: col 1: MASTER NODE (by lag. discretization)
            #         col 2: SLAVE NODE (by lag. discretization)
            #         col 3: length for max. muscle tension
            #         col 4: muscle constant
            #         col 5: hill parameter, a
            #         col 6: hill parameters, b
            #         col 7: force maximum!
    else:
        muscles_info = np.zeros((1,1))  #just to pass placeholder into 
            # "please_Find_Lagrangian_Forces_On_Eulerian_grid function"


    
    # READ IN MUSCLES (IF THERE ARE MUSCLES) #
    if hill_3_muscles_Yes:
        print('  - MUSCLE MODEL (3 Element Hill Model)\n')
        muscles3_info = read_Hill_3Muscle_Points(struct_name)
            #         muscles: col 1: MASTER NODE (by lag. discretization)
            #         col 2: SLAVE NODE (by lag. discretization)
            #         col 3: length for max. muscle tension
            #         col 4: muscle constant
            #         col 5: hill parameter, a
            #         col 6: hill parameters, b
            #         col 7: force maximum!
    else:
        muscles3_info = np.zeros((1,1))  #just to pass placeholder into "please_Find_Lagrangian_Forces_On_Eulerian_grid function"
    
    

    # READ IN USER-DEFINED FORCE MODEL PARAMETERS (IF THERE IS A USER-DEFINED FORCE) #
    if general_force_Yes:
        print('  - GENERAL FORCE MODEL (user-defined force term)\n')
        gen_force_info = read_General_Forcing_Function(struct_name)
        #
        #           
        #   ALL PARAMETERS / FORCE FUNCTION SET BY USER!            
        #            
        #            
    else:
        gen_force_info = 0  

        

    # CONSTRUCT GRAVITY INFORMATION (IF THERE IS GRAVITY) #
    if gravity_Yes:    
        xG = Lag_Struct_Params[14]      # x-Component of Gravity Vector
        yG = Lag_Struct_Params[15]      # y-Component of Gravity Vector
        normG = sqrt( xG**2 + yG**2 )
        gravity_Info = [gravity_Yes, xG/normG, yG/normG]
        #   col 1: flag if considering gravity
        #   col 2: x-component of gravity vector (normalized)
        #   col 3: y-component of gravity vector (normalized)
        
        del xG, yG, normG
        
    else:
        gravity_Info = np.zeros((1,1))


    #
    # BACKGROUND FLOW ITEMS
    #
    print('\n\n--> Background Flow Items\n')
    if ( tracers_Yes == 0 ) and (concentration_Yes == 0):
        print('      (No tracers nor other passive scalars immersed in fluid)\n\n')


    # READ IN TRACERS (IF THERE ARE TRACERS) #
    if tracers_Yes:
        print('  -Tracer Particles included\n')
        nulvar,xT,yT = read_Tracer_Points(struct_name)
        tracers = np.zeros((xT.size,4))
        tracers[0,0] = 1
        tracers[:,1] = xT
        tracers[:,2] = yT
            #tracers_info: col 1: xPt of Tracers
            #              col 2: yPt of Tracers
    else:
        tracers = np.zeros((1,1))



    # READ IN CONCENTRATION (IF THERE IS A BACKGROUND CONCENTRATION) #
    if concentration_Yes:
        print('  -Background concentration included\n')
        C,kDiffusion = read_In_Concentration_Info(struct_name)
            #C:           Initial background concentration
            #kDiffusion:  Diffusion constant for Advection-Diffusion
    else:
        C = 0 # placeholder for plotting 



    # READ IN BRINKMAN TERM (IF THERE IS A BRINKMAN TERM) #
    if brinkman_Yes:
        print('  -Brinkman term included\n')
        Brink,kBrink = read_In_Brinkman_Info(struct_name)
            #Brink:   Initial brinkman permeability matrix
            #kBrink:  Permeability for brinkman model
    else:
        Brink = 0 # placeholder for plotting 


    # CONSTRUCT BOUSSINESQ INFORMATION (IF USING BOUSSINESQ) #
    if boussinesq_Yes:
        print('  -Boussinesq Approximation included\n')
        print('     -> NEED: 1. gravity flag w/ gravity components \n')
        print('              2. background concentration \n')

        if exp_Coeff == 0:
            print('     -> exp_Coeff set to 1.0 by default, was assigned 0 in input2d <-\n')
            exp_Coeff = 1.0

        #if gravity_Info == 0:
        #    print('\n\n\n READ THE ERROR MESSAGE -> YOU MUST FLAG GRAVITY IN INPUT FILE FOR BOUSSINESQ! :) \n\n\n')
        #    error('YOU MUST FLAG GRAVITY IN INPUT FILE FOR BOUSSINESQ! :)')
        #elif concentration_Yes == 0:
        #    print('\n\n\n READ THE ERROR MESSAGE -> YOU MUST HAVE BACKGROUND CONCENTRATION FOR BOUSSINESQ! :) \n\n\n')
        #    error('YOU MUST FLAG CONCENTRATION IN INPUT FILE FOR BOUSSINESQ! :)')

        # Forms Boussinesq forcing terms, e.g., (exp_Coeff)*gVector for Momentum Eq.
        fBouss_X,fBouss_Y = please_Form_Boussinesq_Forcing_Terms(exp_Coeff,Nx,Ny,gravity_Info)

        # Finds initial concentration Laplacian (may not be necessary)
        #Cxx = DD(C,dx,'x')
        #Cyy = DD(C,dy,'y')
        #laplacian_C = Cxx+Cyy
        #C_0 = np.zeros(C.shape) # Define background concentration



    
    # Initialize the initial velocities to zero.
    U = np.zeros((Ny,Nx))                           # x-Eulerian grid velocity
    V = np.zeros((Ny,Nx))                           # y-Eulerian grid velocity
    U.shape
    mVelocity = np.zeros((mass_info.shape[0],2))  # mass-Pt velocity 

    if arb_ext_force_Yes:
        print('  -Artificial External Forcing Onto Fluid Grid\n')
        firstExtForce = 1       # initialize external forcing
        indsExtForce = 0        # initialize for external forcing computation
    
    # ACTUAL TIME-STEPPING IBM SCHEME! 
    #(flags for storing structure connects for printing and printing to .vtk)
    cter = 0 
    ctsave = 0 
    firstPrint = 1 
    loc = 1 
    diffy = 1
    
    # CREATE VIZ_IB2D FOLDER, HIER_IB2D_DATA FOLDER FOR STORING .VTK DATA
    try:
        os.mkdir('viz_IB2d')
    except FileExistsError:
        #File already exists
        pass
    if Output_Params[16] == 1:
        try:
            os.mkdir('hier_IB2d_data')
        except FileExistsError:
            #File already exists
            pass
    os.chdir('viz_IB2d')

    # PRINT BRINKMAN PERMEABILITY TO MATRIX (*if BRINKMAN TERM*) %
    if brinkman_Yes:
        os.chdir('viz_IB2d')
        strNUM = give_String_Number_For_VTK(ctsave)
        brinkName = 'Brink.' + strNUM + '.vtk'
        savevtk_scalar(Brink.T, brinkName, 'Brinkman',dx,dy)
        #clear brinkName strNUM
        os.chdir('..')

    #Save grid_Info and model_Info in human readable format
    #with open('_grid_Info.txt','w') as fobj:
    #    for key in sorted(grid_Info):
    #        fobj.write('{0} = {1}\n'.format(key,grid_Info[key]))
    #with open('_model_Info.txt','w') as fobj:
    #    for key in sorted(model_Info):
    #        fobj.write('{0} = {1}\n'.format(key,model_Info[key]))
    #vizID = open('dumps.visit','w')
    #vizID.write('!NBLOCKS 6\n')
    #os.chdir('..')
    
    #Initialize Vorticity, uMagnitude, and Pressure for initial colormap
    #Print initializations to .vtk
    vort = np.zeros((Nx,Ny)) 
    uMag = np.array(vort) 
    p = np.array(vort)
    lagPts = np.zeros((xLag.size,3))
    lagPts[:,0] = xLag 
    lagPts[:,1] = yLag
    connectsMat,spacing = give_Me_Lag_Pt_Connects(ds,xLag,yLag,Nx,springs_Yes,springs_info)
    Fxh = np.array(vort) 
    Fyh =np.array(vort) 
    F_Lag = np.zeros((xLag.size,2)) 
    print_vtk_files(Output_Params,ctsave,vort,uMag,p,U,V,Lx,Ly,Nx,Ny,lagPts,springs_Yes,\
    connectsMat,tracers,concentration_Yes,C,Fxh,Fyh,F_Lag)
    print('\n |****** Begin IMMERSED BOUNDARY SIMULATION! ******| \n\n')
    print('Current Time(s): {0}\n'.format(current_time))
    ctsave += 1
    
    #
    #
    # * * * * * * * * * * BEGIN TIME-STEPPING! * * * * * * * * * * *
    #
    #
    
    while current_time < T_FINAL:
        
        #
        #
        #******Step 1: Update Position of Boundary of membrane at half time-step ******
        #                 (Variables end with h if it is a half-step)
        #
        xLag_h,yLag_h = please_Move_Lagrangian_Point_Positions(mu, U, V, xLag, yLag,\
            xLag, yLag, x, y, dt/2, grid_Info, 0, poroelastic_Yes, poroelastic_info, F_Poro)
            
        if mass_Yes:
            mass_info, massLagsOld = please_Move_Massive_Boundary(dt/2,\
            mass_info,mVelocity)
           
        if update_Springs_Flag and springs_Yes:
            #This function is application specific, located with main2d
            from update_Springs import update_Springs
            springs_info = update_Springs(dt,current_time,xLag,yLag,springs_info)
            
        if update_Target_Pts and target_pts_Yes:
            #This function is application specific, located with main2d
            from update_Target_Point_Positions import update_Target_Point_Positions
            target_info = update_Target_Point_Positions(dt,current_time,target_info)
            
        if update_Beams_Flag and beams_Yes:
            from update_Beams import update_Beams
            #This function is application specific, located with main2d
            beams_info = update_Beams(dt,current_time,beams_info)

        if update_nonInv_Beams_Flag and nonInv_beams_Yes:
            from update_nonInv_Beams import update_nonInv_Beams
            #This function is application specific, located with main2d
            nonInv_beams_info = update_nonInv_Beams(dt,current_time,beams_info)    

        if update_D_Springs_Flag and d_Springs_Yes:
            from update_Damped_Springs import update_Damped_Springs
            #This function is application specific, located with main2d
            d_springs_info = update_Damped_Springs(dt,current_time,d_springs_info)    
        

        #    
        #
        #*******STEP 2: Calculate Force coming from membrane at half time-step ********
        #
        #
        Fxh, Fyh, F_Mass_Bnd, F_Lag, F_Poro = please_Find_Lagrangian_Forces_On_Eulerian_grid(\
        dt, current_time, xLag_h, yLag_h, xLag_P, yLag_P, x, y, grid_Info, Lag_Struct_Params,\
        springs_info, target_info, beams_info, nonInv_beams_info, muscles_info, muscles3_info,\
        mass_info, d_springs_info, gen_force_info, poroelastic_info)
        
        # Once force is calculated, can finish time-step for massive boundary
        if mass_Yes: #need to test
            # Update Massive Boundary Velocity
            mVelocity_h = please_Update_Massive_Boundary_Velocity(dt/2,mass_info,\
            mVelocity,F_Mass_Bnd,gravity_Info)
            
            # Update Massive Boundary Position for Time-step
            mass_info[:,[1,2]] = massLagsOld
            mass_info,unused = please_Move_Massive_Boundary(dt,mass_info,mVelocity_h)

            # Update Massive Boundary Velocity for Time-step
            mVelocity = please_Update_Massive_Boundary_Velocity(dt,mass_info,\
            mVelocity,F_Mass_Bnd,gravity_Info)
            
        # Add artificial force from fluid boundary, if desired. 
        if arb_ext_force_Yes: #need to test
            # This function is user defined along with main2d
            from please_Compute_External_Forcing import \
            please_Compute_External_Forcing
            # Some of these arguments are mutable. Make sure they are not 
            #   getting assigned to!
            Fx_Arb, Fy_Arb, firstExtForce, indsExtForce = \
            please_Compute_External_Forcing(dt, current_time, x, y, grid_Info,\
            U, V, firstExtForce, indsExtForce)
            Fxh = Fxh + Fx_Arb
            Fyh = Fyh + Fy_Arb


        #  
        #
        #******************* STEP 3: Solve for Fluid motion ************************
        #
        #

        # Add in effect from BOUSSINESQ
        if boussinesq_Yes:
            #Fxh = Fxh + rho*np.dot(fBouss_X,C)
            #Fyh = Fyh + rho*np.dot(fBouss_Y,C)
            Fxh = Fxh + rho*fBouss_X*C #np.dot(fBouss_X,C)
            Fyh = Fyh + rho*fBouss_Y*C #rho*np.dot(fBouss_Y,C)

        # Add in effect from BRINKMAN 
        if brinkman_Yes:
            Fxh = Fxh - mu*Brink*U; # ADD BRINKMAN TERM!
            Fyh = Fyh - mu*Brink*V; # ADD BRINKMAN TERM!
        
        # Update fluid motion
        Uh, Vh, U, V, p =   please_Update_Fluid_Velocity(U, V, Fxh, Fyh, rho, mu, grid_Info, dt, idX, idY)


        
        #
        #
        #**** STEP 4: Update Position of Boundary of membrane again for a half time-step ****
        #
        #
        xLag_P = np.array(xLag_h)   # Stores old Lagrangian x-Values (for muscle model)
        yLag_P = np.array(yLag_h)   # Stores old Lagrangian y-Values (for muscle model)
        #Uh, Vh instead of U,V?
        xLag, yLag = please_Move_Lagrangian_Point_Positions(mu, Uh, Vh, xLag, yLag,\
            xLag_h, yLag_h, x, y, dt, grid_Info,porous_Yes, poroelastic_Yes, poroelastic_info, F_Poro)
            

        # NOTE: SET UP FOR BOTH CLOSED + OPEN SYSTEMS NOW!!!
        if porous_Yes: #need to test
            Por_Mat,nX,nY = please_Compute_Porous_Slip_Velocity(ds,xLag,yLag,\
                porous_info,F_Lag)
            # Update xLag,yLag positions with porous slip velocity
            xLag[porous_info[:,0].astype('int')] = xLag[porous_info[:,0].astype('int')] \
                - dt*Por_Mat[:,0]*nX
            yLag[porous_info[:,0].astype('int')] = yLag[porous_info[:,0].astype('int')] \
                - dt*Por_Mat[:,1]*nY
            # Update storage of xLag, yLag positions inside porous_info structure
            porous_info[:,1] = xLag[porous_info[:,0].astype('int')] #Stores x-Lags as x-Porous Pt. Identities
            porous_info[:,2] = yLag[porous_info[:,0].astype('int')] #Stores y-Lags as y-Porous Pt. Identities    
            xLag = xLag % Lx # If structure goes outside domain
            yLag = yLag % Ly # If structure goes outside domain
            
        # If there are tracers, update tracer positions #
        if tracers_Yes:
            #Uh, Vh instead of U,V?
            xT, yT = please_Move_Lagrangian_Point_Positions(Uh, Vh, xT, yT, xT, yT,\
                x, y, dt, grid_Info,0) #0 for always no porous tracers
            tracers[:,1] = xT
            tracers[:,2] = yT
            
        # If there is a background concentration, update concentration-gradient #
        # Note, C is mutable
        if concentration_Yes: #need to test
           C = please_Update_Adv_Diff_Concentration(C,dt,dx,dy,U,V,kDiffusion)
           
        # Save/Plot Lagrangian/Eulerian Dynamics! #
        if ( ( cter % pDump == 0) and (cter >= pDump) ):
            
            #Compute vorticity, uMagnitude
            vort = give_Me_Vorticity(U,V,dx,dy)
            uMag = give_Me_Magnitude_Velocity(U,V)
            
            #Plot in Matplotlib
            if pMatplotlib:
                loc, diffy = please_Plot_Results(ds,X,Y,U,V,vort,uMag,p,xLag,yLag,\
                    lagPlot,velPlot,vortPlot,pressPlot,uMagPlot,firstPrint,\
                    loc,diffy,spacing)
            
            #Print .vtk files!
            lagPts = np.vstack((xLag, yLag, np.zeros(xLag.size))).T
            print_vtk_files(Output_Params,ctsave,vort,uMag.T,p.T,U.T,V.T,Lx,Ly,Nx,Ny,\
                lagPts,springs_Yes,connectsMat,tracers,concentration_Yes,C,Fxh.T,Fyh.T,F_Lag)
            
            #Print Current Time
            print('Current Time(s): {0:6.6f}\n'.format(current_time))
            
            #Update print counter for filename index
            ctsave+=1
            firstPrint = 0
            
            
        # Update current_time value & counter
        current_time = current_time+dt
        cter += 1
        #wait = input('Press enter to continue...')

    
###########################################################################
#
# FUNCTION: Reads in the # of vertex pts and all the vertex pts from the
#           .vertex file.
#
###########################################################################

def  read_Vertex_Points(struct_name):
    ''' Reads in the number of vertex pts and all the vertex pts from .vertex
    
    Args: 
        struct_name: structure name
        
    Returns:
        N: number of Lagrangian points
        xLag: x-values of Lagrangian mesh
        yLab: y-values of Lagrangian mesh'''

    filename = struct_name+'.vertex'  #Name of file to read in
    with open(filename) as f:
        # First line in the file contains the number of Lagrangian points
        N = int(f.readline().strip())
        # Read in the Lagrangian mesh points
        xLag,yLag = np.loadtxt(f,unpack=True)
       
    return (N,xLag,yLag)
 


 
###########################################################################
#
# FUNCTION: Reads in the # of tracer pts and all the tracer pts from the
#           .tracer file.
#
###########################################################################

def read_Tracer_Points(struct_name):
    '''Reads in the num of tracer pts and all the tracer pts from .tracer
    
    Args:
        struct_name: structure name
        
    Returns:
        N: number of tracer points
        xLag: x-values of tracer points
        yLag: y-values of tracer points'''

    filename = struct_name+'.tracer'  #Name of file to read in
    with open(filename) as f:
        # First line is the number of tracer points
        N = int(f.readline().strip())
        # Now read in the tracer points
        xLag,yLag = np.loadtxt(f,unpack=True)
        
    return (N,xLag,yLag)



    
###########################################################################
#
# FUNCTION: Reads in the diffusion coefficient and initial concentration, C
#
###########################################################################

def read_In_Concentration_Info(struct_name):
    ''' Reads in the diffusion coefficient and initial concentration, C
    
    Args:
        struct_name: structure name
        
    Returns:
        C: concentration
        kDiff: coefficient of diffusion'''

    filename = struct_name+'.concentration'  #Name of file to read in
    
    with open(filename) as f:
        kDiff = float(f.readline().strip()) #first line contains coeff of diff
        C = np.loadtxt(f)

    return (C,kDiff)


###########################################################################
#
# FUNCTION: Reads in the Brinkman permeability, k and initial Brinkman term, Brink
#
###########################################################################

def read_In_Brinkman_Info(struct_name):
    ''' Reads in the diffusion coefficient and initial concentration, C
    
    Args:
        struct_name: structure name
        
    Returns:
        Brink:  initial Brinkman region (1/kBrink)
        kBrink: coefficient of permeability'''

    filename = struct_name+'.brinkman'  #Name of file to read in
    
    with open(filename) as f:
        kBrink = float(f.readline().strip()) #first line contains coeff of diff
        C = np.loadtxt(f)

    return (Brink,kBrink)


    
###########################################################################
#
# FUNCTION: Reads in the # of springs and all MASTER NODEs, SLAVE NODEs,
#           spring STIFFNESSES, spring RESTING LENGTHS
#
###########################################################################

def read_Spring_Points(struct_name):
    ''' Reads in the num of springs, master/slave nodes, spring stiffness/resting lengths, deg. of non-linearity
    
    Args:
        struct_name: structure name
        
    Returns:
        springs: above info stored in columns'''

    filename = struct_name+'.spring'  #Name of file to read in
    try:
        springs = np.genfromtxt(filename,skip_header=1,
            missing_values=['-NaN', '-nan', 'N/A', 'NA', 'NULL', 'NaN', 'nan'],filling_values=1)
    except ValueError:
        print('Failed to load spring data from {}.\n'.format(filename)+
              'Check that all rows (after the first) have the same number of columns\n'+
              'N/A, NA, NULL, NaN, and nan can be used to denote missing values if\n'+
              'linear and non-linear springs are mixed (these entries will be replaced with a 1).')
        raise

    # If no specified degreee on non-linearity in .spring file
    n,m = springs.shape
    if (m<5):
        springs2 = np.ones([n,1])                           # DEFAULT DEG. NL if nothing inputted

        springs = np.concatenate( (springs,springs2), axis=1 )  # Combine springs w/ DEG. NONLINEARITY if not there before

    #springs: col 1: starting spring pt (by lag. discretization)
    #         col 2: ending spring pt. (by lag. discretization)
    #         col 3: spring stiffness
    #         col 4: spring resting lengths
    #         col 5: degree of non-linearity
    
    return springs



################################################################################
#
# FUNCTION: Reads in the # of DAMPED springs and all MASTER NODEs, SLAVE NODEs,
#           spring STIFFNESSES, spring RESTING LENGTHS, spring DAMPING coeffs.
#
################################################################################

def read_Damped_Spring_Points(struct_name):
    ''' Reads in the num of springs, master/slave nodes, spring stiffness/resting lengths
    
    Args:
        struct_name: structure name
        
    Returns:
        springs: above info stored in columns'''

    filename = struct_name+'.d_spring'  #Name of file to read in
    with open(filename) as f:
    #Store elements on .spring file into a matrix starting w/ 2nd row of read in data.
        dSprings = np.loadtxt(f,skiprows=1,usecols=(0,1,2,3,4))

    #springs: col 1: starting spring pt (by lag. discretization)
    #         col 2: ending spring pt. (by lag. discretization)
    #         col 3: spring stiffness
    #         col 4: spring resting lengths
    #         col 5: damping coefficient
    
    return dSprings

    
###########################################################################
#
# FUNCTION: Reads in the # of muscles and all MASTER NODEs, SLAVE NODEs,
#           length for max. muscle tension, muscle constant, hill
#           parameters (a and b), and Force-Max
#
###########################################################################

def read_Muscle_Points(struct_name):
    ''' Reads in the # of muscles and other information
    
    Args:
        struct_name: structure name
        
    Returns:
        muscles: array of muscle info'''

    filename = struct_name+'.muscle'  #Name of file to read in
    with open(filename) as f:
        muscles = np.loadtxt(f,skiprows=1,usecols=(0,1,2,3,4,5,6))

    #muscles: col 1: MASTER NODE (by lag. discretization)
    #         col 2: SLAVE NODE (by lag. discretization)
    #         col 3: length for max. muscle tension
    #         col 4: muscle constant
    #         col 5: hill parameter, a
    #         col 6: hill parameters, b
    #         col 7: force maximum!
    
    return muscles
    



###########################################################################
#
# FUNCTION: Reads in the # of muscles and all MASTER NODEs, SLAVE NODEs,
#           length for max. muscle tension, muscle constant, hill
#           parameters (a and b), and Force-Max
#
###########################################################################

def read_Hill_3Muscle_Points(struct_name): #untested
    ''' Reads in the num of muscles and other info
    
    Args:
        struct_name: structure name
        
    Returns:
        muscles: muscle info'''

    filename = struct_name+'.muscle'  #Name of file to read in
    with open(filename) as f:
        muscles = np.loadtxt(f,skiprows=1,usecols=(0,1,2,3,4,5,6))

    #muscles: col 1: MASTER NODE (by lag. discretization)
    #         col 2: SLAVE NODE (by lag. discretization)
    #         col 3: length for max. muscle tension
    #         col 4: muscle constant
    #         col 5: hill parameter, a
    #         col 6: hill parameters, b
    #         col 7: force maximum!
    
    return muscles




    
###########################################################################
#
# FUNCTION: Reads in the # of MASS PTS, Mass-STIFFNESSES, and Mass-VALUE
#           
#
###########################################################################

def read_Mass_Points(struct_name):
    ''' Reads in the num of mass pts, mass-spring stiffness, and mass-value
    
    Args:
        struct_name: structure name
        
    Returns:
        masses: array of mass info'''

    filename = struct_name+'.mass'  #Name of file to read in
    with open(filename) as f:
        masses = np.loadtxt(f,skiprows=1,usecols=(0,1,2))

    #masses:  col 1: Lag Pt. ID w/ Associated Mass Pt.
    #         col 2: "Mass-Spring" stiffness Parameter
    #         col 3: Mass Value Parameter
    
    return masses





###########################################################################
#
# FUNCTION: Reads in the # of TARGET PTS, TARGET-PT-NODEs, and their
#           Target-STIFFNESSES
#
###########################################################################

def read_Target_Points(struct_name):
    ''' Reads in the num of target pts, target-pt-nodes, and target-stiffness
    
    Args:
        struct_name: structure name
        
    Returns:
        targets: array of target info'''
    filename = struct_name+'.target'  #Name of file to read in
    with open(filename) as f:
        targets = np.loadtxt(f,skiprows=1,usecols=(0,1))

    #targets: col 1: Lag Pt. ID w/ Associated Target Pt.
    #         col 2: target STIFFNESSES
    
    return targets
    

    
###########################################################################
#
# FUNCTION: Reads in the # of POROUS PTS, POROUS-PT-NODEs, and their
#           POROUSITY-COEFFICIENTS
#
###########################################################################

def read_Porous_Points(struct_name):
    ''' Reads in the num of porous pts, pt-nodes, and porousity coefficients
    
    Args:
        struct_name: structure name
        
    Returns:
        porosity: array of porosity info'''

    filename = struct_name+'.porous'  #Name of file to read in
    with open(filename) as f:
        porosity = np.loadtxt(f,skiprows=1,usecols=(0,1,2))

    #porous:  col 1: Lag Pt. ID w/ Associated Porous Pt.
    #         col 2: Porosity coefficient
    
    return porosity



###########################################################################
#
# FUNCTION: Reads in the # of POROUS PTS, POROUS-PT-NODEs, and their
#           POROUSITY-COEFFICIENTS
#
###########################################################################

def read_PoroElastic_Points(struct_name):
    ''' Reads in the num of porous pts, pt-nodes, and porousity coefficients
    
    Args:
        struct_name: structure name
        
    Returns:
        porosity: array of porosity info'''

    filename = struct_name+'.poroelastic'  #Name of file to read in
    with open(filename) as f:
        porosity = np.loadtxt(f,skiprows=1,usecols=(0,1))

    #porous:  col 1: Lag Pt. ID w/ Associated Porous Pt.
    #         col 2: Brinkman coefficient
    
    return porosity










###########################################################################
#
# FUNCTION: Reads in the # of beams and all 1st Pt, MIDDLE Pt, and 3rd Pt
#           beam STIFFNESSES, and CURVATURE
#
###########################################################################

def read_Beam_Points(struct_name):
    ''' Reads in the num of beams, 1st pt, middle pt, 3 pt stiffness, curvature
    
    Args:
        struct_name: structure name
        
    Returns:
        beams: array of beam info'''

    filename = struct_name+'.beam'  #Name of file to read in
    with open(filename) as f:
        beams = np.loadtxt(f,skiprows=1,usecols=(0,1,2,3,4))

    #beams:      col 1: 1ST PT.
    #            col 2: MIDDLE PT. (where force is exerted)
    #            col 3: 3RD PT.
    #            col 4: beam stiffness
    #            col 5: curavture
    
    return beams


###########################################################################
#
# FUNCTION: Reads in the # of beams and all 1st Pt, MIDDLE Pt, and 3rd Pt
#           beam STIFFNESSES, and CURVATURE for NON-INVARIANT BEAMS
#
###########################################################################

def read_nonInv_Beam_Points(struct_name):
    ''' Reads in the num of beams, 1st pt, middle pt, 3 pt stiffness, curvature
    
    Args:
        struct_name: structure name
        
    Returns:
        beams: array of beam info'''

    filename = struct_name+'.nonInv_beam'  #Name of file to read in
    with open(filename) as f:
        beams = np.loadtxt(f,skiprows=1,usecols=(0,1,2,3,4,5))

    #beams:      col 1: 1ST PT.
    #            col 2: MIDDLE PT. (where force is exerted)
    #            col 3: 3RD PT.
    #            col 4: beam stiffness
    #            col 5: curavture
    
    return beams   


###########################################################################
#
# FUNCTION: READS IN ALL THE DATA FOR THE USER-DEFINED FORCE FUNCTION!!!
#           NOTE: DOES NOT SPECIFY HOW MANY PARAMETERS THERE ARE.
#           NOTE: COMPLETELY USER-DEFINED
#
###########################################################################

def  read_General_Forcing_Function(struct_name):
    ''' Reads in all data from user defined force .user_force file
    
    Args: 
        struct_name: structure name
        
    Returns:
        All parameters set by user for user-defined force'''

    filename = struct_name+'.user_force'  #Name of file to read in
    with open(filename) as f:
        # First line in the file contains the number of Lagrangian points
        N = int(f.readline().strip())
        # Read in the Lagrangian mesh points
        force_general = np.loadtxt(f,unpack=True)

    return force_general   
    
    
##############################################################################
#
# FUNCTION: give me Connects Vector for printing Lagrangian .vtk info!
#
##############################################################################

def give_Me_Lag_Pt_Connects(ds,xLag,yLag,Nx,springs_Yes,springs_info):
    ''' Give me Connects Vector for printing Lagrangian .vtk info!

    Args:
        ds:
        xLag:
        yLag:
        Nx:
        
    Returns:
        connectsMat:
        space:'''

    #springs_info: col 1: starting spring pt (by lag. discretization)
    #              col 2: ending spring pt. (by lag. discretization)
    #              col 3: spring stiffness
    #              col 4: spring resting lengths

    space = 20*ds

    nSprings = springs_info.shape
    connectsMat = np.zeros([nSprings[0],2])

    if springs_Yes:
        connectsMat[:,0] = springs_info[:,0] # (for .vtk counting)
        connectsMat[:,1] = springs_info[:,1] # (for .vtk counting)
    else:
        connectsMat = 0

    #N = xLag.size

    #if Nx <= 32:
    #    space = 5*ds
    #elif Nx <= 64:
    #   space = 5*ds
    #elif Nx <=128:
    #   space = 5*ds
    #elif Nx <=256:
    #    space = 10*ds
    #elif Nx <= 512:
    #    space = 20*ds
    #else:
    #    space = 40*ds
        
    #dist = np.zeros(N)
    #dist[:-1] = np.sqrt( (xLag[:-1]-xLag[1:])**2 + (yLag[:-1]-yLag[1:])**2 )
    #dist[-1] = sqrt( (xLag[-1]-xLag[0])**2 + (yLag[-1]-yLag[0])**2 )
    #collect the indices where dist < space, for Cpp notation (and .vtk counting)
    #connectsMat0 = np.where(dist<space)[0] #always returns tuple of arrays of
                                           # indices, so need first array in tuple
    #if connectsMat0.size > 0:
    #    connectsMat1 = (connectsMat0 + 1) % N # N index should wrap back to 0
    #connectsMat = np.vstack((connectsMat0,connectsMat1)).T
    
    return (connectsMat,space)


##############################################################################
#
# FUNCTION: gives appropriate string number for filename in printing the
# .vtk files.
#
##############################################################################

def print_vtk_files(Output_Params,ctsave,vort,uMag,p,U,V,Lx,Ly,Nx,Ny,lagPts,springs_Yes,\
    connectsMat,tracers,concentration_Yes,C,fXGrid,fYGrid,F_Lag):
    ''' Gives appropriate string number for filename in printing the .vtk files'''


    #               Output_Params[0]: print_dump
    #                            [1]: plot_Matlab
    #                            [2]: plot_LagPts
    #                            [3]: plot_Velocity
    #                            [4]: plot_Vorticity
    #                            [5]: plot_MagVelocity
    #                            [6]: plot_Pressure
    #                            [7]:  save_Vorticity 
    #                            [8]:  save_Pressure 
    #                            [9]: save_uVec 
    #                            [10]: save_uMag 
    #                            [11]: save_uX 
    #                            [12]: save_uY 
    #                            [13]: save_fMag 
    #                            [14]: save_fX 
    #                            [15]: save_fY 
    #                            [16]: save_hier 


    #Give EULERIAN spacing for grid
    dx = Lx/Nx 
    dy = Ly/Ny


    #Go into viz_IB2d directory. This was throwing an error because we're already there!
    if os.path.split(os.getcwd())[1] != 'viz_IB2d':
        os.chdir('viz_IB2d')

    #Find string number for storing files
    strNUM = give_String_Number_For_VTK(ctsave)
    lagPtsName = 'lagsPts.'+strNUM+'.vtk'

    #Print Lagrangian Pts to .vtk format
    savevtk_points(lagPts, lagPtsName, 'lagPts')

    # Print Spring Connections (* if springs *)
    if springs_Yes:
        #Print Lagrangian Pts w/ CONNECTIONS to .vtk format
        lagPtsConName = 'lagPtsConnect.'+strNUM+'.vtk'
        savevtk_points_connects(lagPts, lagPtsConName, 'lagPtsConnected',connectsMat )

    #Print Tracer Pts (*if tracers*)
    if tracers[0,0] == 1:
        tracersPtsName = 'tracer.'+strNUM+'.vtk'
        #tMatrix = tracers[:,1:4]
        savevtk_points(tracers[:,1:4],tracersPtsName, 'tracers') 


    #Print SCALAR DATA (i.e., colormap data) to .vtk file
    if Output_Params[7] == 1:
        vortfName = 'Omega.'+strNUM+'.vtk'
        savevtk_scalar(vort, vortfName, 'Omega',dx,dy)

    if Output_Params[8] == 1:
        pfName = 'P.'+strNUM+'.vtk'
        savevtk_scalar(p, pfName, 'P',dx,dy)

    if Output_Params[10] == 1:
        uMagfName = 'uMag.'+strNUM+'.vtk'
        savevtk_scalar(uMag, uMagfName, 'uMag',dx,dy)

    if Output_Params[11] == 1:
        uXName = 'uX.'+strNUM+'.vtk'
        savevtk_scalar(U, uXName, 'uX',dx,dy)
    
    if Output_Params[12] == 1:
        uYName = 'uY.'+strNUM+'.vtk'    
        savevtk_scalar(V, uYName, 'uY',dx,dy)

    if Output_Params[14] == 1:
        fXName = 'fX.'+strNUM+'.vtk'
        savevtk_scalar(fXGrid, fXName, 'fX',dx,dy)

    if Output_Params[15] == 1:
        fYName = 'fY.'+strNUM+'.vtk'
        savevtk_scalar(fYGrid, fYName, 'fY',dx,dy)
    
    if Output_Params[13] == 1:
        fMagName = 'fMag.'+strNUM+'.vtk'
        savevtk_scalar(np.sqrt( fXGrid*fXGrid + fYGrid*fYGrid ), fMagName, 'fMag',dx,dy)

    if concentration_Yes:
        confName = 'concentration.'+strNUM+'.vtk'
        savevtk_scalar(C.T, confName, 'Concentration',dx,dy)

    #Print VECTOR DATA (i.e., velocity data) to .vtk file
    if Output_Params[9] == 1:
        velocityName = 'u.'+strNUM+'.vtk'
        savevtk_vector(U, V, velocityName, 'u',dx,dy)

    #Get out of viz_IB2d folder
    os.chdir('..')


    #
    # Print Lagrangian Force Data to hier_IB2d_data folder
    #
    NLagPts = lagPts.shape[1]
    
    if Output_Params[16] == 1:

    # THE CASE IF LESS THAN (or =) to 5 Lag. Pts. 
        if NLagPts <= 5:
            os.chdir('hier_IB2d_data') #change directory to hier-data folder
        
            # Save x-y force data!
            fLag_XName = 'fX_Lag.'+strNUM+'.vtk'
            fLag_YName = 'fY_Lag.'+strNUM+'.vtk'
            savevtk_points_with_scalar_data( lagPts, F_Lag[:,0], fLag_XName, 'fX_Lag')
            savevtk_points_with_scalar_data( lagPts, F_Lag[:,1], fLag_YName, 'fY_Lag')

            # Define force magnitude name
            fMagName = 'fMag.'+strNUM+'.vtk'
        
            # Compute magnitude of forces on Lagrangian boundary
            fLagMag = np.sqrt( F_Lag[:,0]*F_Lag[:,0] + F_Lag[:,1]*F_Lag[:,1] ) 
        
            # Print UNSTRUCTURED POINT DATA w/ SCALAR associated with it
            savevtk_points_with_scalar_data( lagPts, fLagMag, fMagName, 'fMag')
        
            # Get out of hier_IB2d_data folder
            os.chdir('..') 
    
        # THE CASE IF GREATER THAN 5 Lag. Pts. 
        else:
            F_Tan_Mag,F_Normal_Mag = please_Compute_Normal_Tangential_Forces_On_Lag_Pts(lagPts,F_Lag)

            os.chdir('hier_IB2d_data') #change directory to hier-data folder

            # Save x-y force data!
            fLag_XName = 'fX_Lag.'+strNUM+'.vtk'
            fLag_YName = 'fY_Lag.'+strNUM+'.vtk'
            savevtk_points_with_scalar_data( lagPts, F_Lag[:,0], fLag_XName, 'fX_Lag')
            savevtk_points_with_scalar_data( lagPts, F_Lag[:,1], fLag_YName, 'fY_Lag')

            # Save force magnitude, mag. normal force, and mag. tangential force
            fMagName = 'fMag.'+strNUM+'.vtk'
            fNormalName = 'fNorm.'+strNUM+'.vtk'
            fTangentName = 'fTan.'+strNUM+'.vtk'

            # Compute magnitude of forces on Lagrangian boundary
            fLagMag = np.sqrt( F_Lag[:,0]*F_Lag[:,0] + F_Lag[:,1]*F_Lag[:,1] )

            # Print UNSTRUCTURED POINT DATA w/ SCALAR associated with it
            savevtk_points_with_scalar_data( lagPts, fLagMag, fMagName, 'fMag')
            savevtk_points_with_scalar_data( lagPts, F_Normal_Mag, fNormalName, 'fNorm')
            savevtk_points_with_scalar_data( lagPts, F_Tan_Mag, fTangentName, 'fTan')

            # Get out of hier_IB2d_data folder
            os.chdir('..') 
    
    
##############################################################################
#
# FUNCTION: gives appropriate string number for filename in printing the
# .vtk files.
#
##############################################################################

def give_String_Number_For_VTK(num):
    ''' Gives appropriate string number for filename in printing .vtk files
    
    Args:
        num: number of file to be printed
        
    Returns:
        strNUM: string number for filename'''

    if num < 10:
        strNUM = '000'+str(num)
    elif num < 100:
        strNUM = '00'+str(num)
    elif num<1000:
        strNUM = '0'+str(num)
    else:
        strNUM = str(num)
        
    return strNUM

##############################################################################
#
# FUNCTION: prints matrix vector data to vtk formated file
#
##############################################################################
def savevtk_points_connects( X, filename, vectorName,connectsMat):
    '''Prints matrix vector data to vtk formated file
    
    Args:
        X: Matrix of size Nx3
        filename: File name
        vectorname:
        connectsMat:'''

    N = X.shape[0]
    Nc = connectsMat.shape[0]

    if C_flag:
        #Just add the measure of time for transforming the 
        nX = np.ascontiguousarray(X, dtype=np.float64)
        nconnectsMat = np.ascontiguousarray(connectsMat, dtype=np.float64)
        write.savevtk_points_connects_write(N,Nc,nX,filename,vectorName,nconnectsMat)
    elif vtk_lib_flag:
        ### UNSTRUCTURED_GRID - POINTS, CELLS, CELL_TYPES ###
        # Create Points data type
        vtk_points = vtk.vtkPoints()
        vtk_points.SetData(numpy_support.numpy_to_vtk(
                           np.require(X, dtype=np.float64, requirements=['C']),
                           deep=1))
        # Create CellArray data type
        vtk_cells = vtk.vtkCellArray()
        # Create connectivity list
        cells = np.ones((Nc,3), dtype=np.int64)
        cells[:,0] *= connectsMat.shape[1]
        cells[:,1:] = connectsMat.astype(np.int64)
        cells = np.ravel(cells)
        # Create cell structure with args: number of cells, connectivity list
        vtk_cells.SetCells(Nc, numpy_support.numpy_to_vtkIdTypeArray(cells))
        vtk_obj = vtk.vtkUnstructuredGrid()
        # Set cells with cell types, cell structure
        vtk_obj.SetCells(list(np.ones(Nc, dtype=np.int64)*3), vtk_cells)
        vtk_obj.SetPoints(vtk_points)
        # write out
        writer = vtk.vtkGenericDataObjectWriter()
        writer.SetFileName(filename)
        writer.SetHeader(vectorName+" generated by IB2d Python")
        writer.SetInputDataObject(vtk_obj)
        writer.Update()
        writer.Write()
    else:
        with open(filename,'w') as file:
            file.write('# vtk DataFile Version 2.0\n')
            file.write(vectorName+'\n')
            file.write('ASCII\n')
            file.write('DATASET UNSTRUCTURED_GRID\n\n')
            #
            file.write('POINTS {0} float\n'.format(N))
            for ii in range(N):
                file.write('{0:.15e} {1:.15e} {2:.15e}\n'.format(X[ii,0],X[ii,1],X[ii,2]))
            file.write('\n')
            #
            #First: # of "Cells", Second: Total # of info inputed following
            file.write('CELLS {0} {1}\n'.format(Nc,3*Nc))
            for s in range(Nc):
                file.write('{0} {1:d} {2:d}\n'.format(2,int(connectsMat[s,0]),int(connectsMat[s,1]) ) ) 
            file.write('\n')
            #
            file.write('CELL_TYPES {0}\n'.format(Nc)) # N = # of "Cells"
            for ii in range(Nc):
               file.write('3 ')
            file.write('\n')
    



##############################################################################
#
# FUNCTION: prints matrix vector data to vtk formated file
#
##############################################################################
def savevtk_points( X, filename, vectorName):
    ''' Prints matrix vector data to vtk formated file
    
    Args:
        X: Matrix of size Nx3
        filename:
        vectorName:'''

    N = X.shape[0]

    if C_flag:
        nX = np.ascontiguousarray(X, dtype=np.float64)
        write.savevtk_points_write(N,nX,filename,vectorName)
    elif vtk_lib_flag:
        ### UNSTRUCTURED_GRID - POINTS, CELLS, CELL_TYPES ###
        # Create Points data type
        vtk_points = vtk.vtkPoints()
        vtk_points.SetData(numpy_support.numpy_to_vtk(
                           np.require(X, dtype=np.float64, requirements=['C']),
                           deep=1))
        # Create CellArray data type
        vtk_cells = vtk.vtkCellArray()
        # Create connectivity list
        cells = np.stack((np.ones(N,dtype=np.int64), np.arange(N,dtype=np.int64)),
                         axis=-1)
        cells_r = np.ravel(cells)
        # Create cell structure with args: number of cells, connectivity list
        vtk_cells.SetCells(N, numpy_support.numpy_to_vtkIdTypeArray(cells_r))
        vtk_obj = vtk.vtkUnstructuredGrid()
        # Set cells with cell types, cell structure
        if cells.shape[0] > 1:
            vtk_obj.SetCells(list(cells[:,0]), vtk_cells)
        else:
            vtk_obj.SetCells(cells[0,0], vtk_cells)
        vtk_obj.SetPoints(vtk_points)
        # write out
        writer = vtk.vtkGenericDataObjectWriter()
        writer.SetFileName(filename)
        writer.SetHeader(vectorName+" generated by IB2d Python")
        writer.SetInputDataObject(vtk_obj)
        writer.Update()
        writer.Write()
    else:
        with open(filename,'w') as file:
            file.write('# vtk DataFile Version 2.0\n')
            file.write(vectorName+'\n')
            file.write('ASCII\n')
            file.write('DATASET UNSTRUCTURED_GRID\n\n')
            file.write('POINTS {0} float\n'.format(N))
            for ii in range(N):
                file.write('{0:.15e} {1:.15e} {2:.15e}\n'.format(X[ii,0],X[ii,1],X[ii,2]))
            file.write('\n')
            #
            #First: # of "Cells", Second: Total # of info inputed following
            file.write('CELLS {0} {1}\n'.format(N,2*N))
            for s in range(N):
                file.write('{0} {1}\n'.format(1,s))
            file.write('\n')
            #
            file.write('CELL_TYPES {0}\n'.format(N)) # N = # of "Cells"
            for ii in range(N):
               file.write('1 ')
            file.write('\n')



##############################################################################
#
# FUNCTION: prints matrix vector data to vtk formated file
#
##############################################################################
def savevtk_vector(X, Y, filename, vectorName,dx,dy):
    ''' Prints matrix vector data to vtk formated file.
    
    Args:
        X: 2-D ndarray
        Y: 2-D ndarray
        filename: file name
        vectorName: name of the variable
        dx:
        dy:'''
    #  ?? Legacy:
    #  savevtkvector Save a 3-D vector array in VTK format
    #  savevtkvector(X,Y,Z,filename) saves a 3-D vector of any size to
    #  filename in VTK format. X, Y and Z should be arrays of the same
    #  size, each storing speeds in the a single Cartesian directions.
    
    #  Christopher's note:
    #   3-D is clearly broken in this code, but there were still some reminants 
    #   in the matlab version. Given the choice of doing try/except blocks to
    #   keep these reminants or to kill them entirely, I'm choosing to kill them.
    #   So, specifically, nz is now gone.
    
    assert (X.shape == Y.shape), 'Error: velocity arrays of unequal size'
    nx, ny = X.shape
    

    XRow = X.shape[0]
    XCol = X.shape[1]
    YRow = Y.shape[0]
    YCol = Y.shape[1]


    if C_flag:
        nX = np.ascontiguousarray(X, dtype=np.float64)
        nY = np.ascontiguousarray(Y, dtype=np.float64)
        write.savevtk_vector(XRow,XCol,YRow,YCol,nX,nY,filename,vectorName,dx,dy)
    elif vtk_lib_flag:
        ### STRUCTURED_POINTS - VECTORS ###
        # Collect info to write
        origin = (0.0, 0.0, 0.0)
        spacing = (dx, dy, 1.0)
        dimensions = (nx, ny, 1)
        vec_array = np.require(np.stack((np.ravel(X,order='F'),
                               np.ravel(Y,order='F'), np.zeros(Y.size)), axis=-1),
                               requirements=['C'])
        vtk_array = numpy_support.numpy_to_vtk(vec_array)
        vtk_array.SetName(vectorName)
        # Create data object
        vtk_obj = vtk.vtkStructuredPoints()
        vtk_obj.SetOrigin(origin)
        vtk_obj.SetSpacing(spacing)
        vtk_obj.SetDimensions(dimensions)
        vtk_obj.GetPointData().SetVectors(vtk_array)
        # Write out data
        writer = vtk.vtkGenericDataObjectWriter()
        writer.SetFileName(filename)
        writer.SetHeader("Generated by IB2d Python")
        writer.SetInputDataObject(vtk_obj)
        writer.Update()
        writer.Write()
    else:
        with open(filename,'w') as fid:
            fid.write('# vtk DataFile Version 2.0\n')
            fid.write('Comment goes here\n')
            fid.write('ASCII\n')
            fid.write('\n')
            fid.write('DATASET STRUCTURED_POINTS\n')
            # 1 below was nz
            fid.write('DIMENSIONS    {0}   {1}   {2}\n'.format(ny, nx, 1))
            fid.write('\n')
            fid.write('ORIGIN    0.000   0.000   0.000\n')
            #fid.write('SPACING   1.000   1.000   1.000\n') #if want [1,32]x[1,32] rather than [0,Lx]x[0,Ly]
            fid.write('SPACING   '+str(dx)+str(' ')+str(dy)+'   1.000\n')
            fid.write('\n')
            fid.write('POINT_DATA   {0}\n'.format(nx*ny))
            fid.write('VECTORS '+vectorName+' double\n')
            fid.write('\n')
            for b in range(ny):
                for c in range(nx):
                    fid.write('{0} '.format(X[c,b]))
                    fid.write('{0} '.format(Y[c,b]))
                    fid.write('0 ')
                fid.write('\n')
    #Python 3.5 automatically opens in text mode unless otherwise specified


##############################################################################
#
# FUNCTION: prints scalar matrix to vtk formated file
#
##############################################################################
def savevtk_scalar(array, filename, dataName,dx,dy):
    ''' Prints scalar matrix to vtk formatted file.
    
    Args:
        array: 2-D ndarray
        filename: file name
        dataName: string describing the data
        dx:
        dy:'''

    #  ?? Legacy:
    #  savevtk Save a 3-D scalar array in VTK format.
    #  savevtk(array, filename) saves a 3-D array of any size to
    #  filename in VTK format.
    
    #  Christopher's note:
    #   3-D is clearly broken in this code, but there were still some reminants 
    #   in the matlab version. Given the choice of doing try/except blocks to
    #   keep these reminants or to kill them entirely, I'm choosing to kill them.
    #   So, specifically, nz is now gone.
    ny,nx = array.shape
    if C_flag:
        narray = np.ascontiguousarray(array, dtype=np.float64)
        write.savevtk_scalar(ny,nx,narray,filename,dataName,dx,dy)
    elif vtk_lib_flag:
        ### STRUCTURED_POINTS - SCALARS ###
        # Collect info to write
        origin = (0.0, 0.0, 0.0)
        spacing = (dx, dy, 1.0)
        dimensions = (ny, nx, 1)
        array_r = np.ravel(array,order='F')
        vtk_array = numpy_support.numpy_to_vtk(array_r)
        vtk_array.SetName(dataName)
        # Create data object
        vtk_obj = vtk.vtkStructuredPoints()
        vtk_obj.SetOrigin(origin)
        vtk_obj.SetSpacing(spacing)
        vtk_obj.SetDimensions(dimensions)
        vtk_obj.GetPointData().SetScalars(vtk_array)
        # Write out data
        writer = vtk.vtkGenericDataObjectWriter()
        writer.SetFileName(filename)
        writer.SetHeader("Generated by IB2d Python")
        writer.SetInputDataObject(vtk_obj)
        writer.Update()
        writer.Write()
    else:
        with open(filename,'w') as fid:
            fid.write('# vtk DataFile Version 2.0\n')
            fid.write('Comment goes here\n')
            fid.write('ASCII\n')
            fid.write('\n')
            fid.write('DATASET STRUCTURED_POINTS\n')
            # 1 below was nz
            fid.write('DIMENSIONS    {0}   {1}   {2}\n'.format(ny, nx, 1))
            fid.write('\n')
            fid.write('ORIGIN    0.000   0.000   0.000\n')
            #fid.write('SPACING   1.000   1.000   1.000\n') #if want [1,32]x[1,32] rather than [0,Lx]x[0,Ly]
            fid.write('SPACING   '+str(dx)+str(' ')+str(dy)+'   1.000\n')
            fid.write('\n')
            # The 1 below was nz
            fid.write('POINT_DATA   {0}\n'.format(nx*ny*1))
            fid.write('SCALARS '+dataName+' double\n')
            fid.write('LOOKUP_TABLE default\n')
            fid.write('\n')
            for b in range(nx):
                for c in range(ny):
                    fid.write('{0} '.format(array[c,b]))
                fid.write('\n')
        #Python 3.5 automatically opens in text mode unless otherwise specified
            

##############################################################################
#
# FUNCTION: prints unstructured point data w/ associated scalar value to vtk 
#           formated file
#
##############################################################################

def savevtk_points_with_scalar_data( X, scalarArray, filename, dataName):
    ''' Prints matrix vector data to vtk formated file
    
    Args:
        X: Matrix of size Nx3
        scalarArray:
        filename:
        dataName: string describing the data'''

    # X is matrix of size Nx3 
    #              Col 1: x-data
    #              Col 2: y-data
    #              Col 3: z-data
    # scalarArray: Scalar array you are assigning to each point
    # filename:    What you are saving the VTK file as (string)
    # dataName:  What you are naming the data you're printing (string)

    N = X.shape[0]
    nx = scalarArray.shape[0]


    if C_flag:
        nX = np.ascontiguousarray(X, dtype=np.float64)
        write.savevtk_points_write(N,nX,filename,dataName)
    elif vtk_lib_flag:
        ### UNSTRUCTURED_GRID - POINTS, POINT_DATA (SCALARS) ###
        # Create Points data type
        vtk_points = vtk.vtkPoints()
        vtk_points.SetData(numpy_support.numpy_to_vtk(
                           np.require(X, dtype=np.float64, requirements=['C']),
                           deep=1))
        scalarArray_r = np.require(scalarArray, dtype=np.float64, requirements=['C'])
        vtk_scalarArray = numpy_support.numpy_to_vtk(scalarArray_r)
        vtk_scalarArray.SetName(dataName)
        vtk_obj = vtk.vtkUnstructuredGrid()
        vtk_obj.SetPoints(vtk_points)
        vtk_obj.GetPointData().SetScalars(vtk_scalarArray)
        # write out
        writer = vtk.vtkGenericDataObjectWriter()
        writer.SetFileName(filename)
        writer.SetHeader(dataName+" generated by IB2d Python")
        writer.SetInputDataObject(vtk_obj)
        writer.Update()
        writer.Write()
    else:
        with open(filename,'w') as file:
            file.write('# vtk DataFile Version 2.0\n')
            file.write(dataName+'\n')
            file.write('ASCII\n')
            file.write('DATASET UNSTRUCTURED_GRID\n\n')
            file.write('POINTS {0} float\n'.format(N))
            for ii in range(N):
                file.write('{0:.15e} {1:.15e} {2:.15e}\n'.format(X[ii,0],X[ii,1],X[ii,2]))
            file.write('\n')
            #
            file.write('POINT_DATA   {0}\n'.format(nx*1*1))
            file.write('SCALARS '+dataName+' double\n')
            file.write('LOOKUP_TABLE default\n')
            file.write('\n')
            for c in range(nx):
                file.write('{0} '.format(scalarArray[c]))
                file.write('\n')
            #Python 3.5 automatically opens in text mode unless otherwise specified


            
##############################################################################
#
# FUNCTION: Computes vorticity from two matrices, U and V, where each
# matrix is the discretized field of velocity values either for x or y,
# respectively.
#
##############################################################################

def give_Me_Vorticity(U,V,dx,dy):
    ''' Computes vorticity from two matrices, discretized field of velocity values
    
    Args:
        U: x-directed vorticity
        V: y-directed vorticity
        dx:
        dy:
        
    Returns:
        vort:'''

    # w = ( dv/dx - du/dy )\hat{z}

    #Compute dv/dx using central differencing! (maintains periodicity)
    dvdx = D(V,dx,'x')

    #Compute du/dy using central differencing! (maintains periodicity)
    dudy = D(U,dy,'y')

    #Compute vorticity
    vort = ( dvdx - dudy )

    #Take transpose so all lines up
    vort = vort.T

    return vort


##############################################################################
#
# FUNCTION: Computes velocity magnitude from two matrices, U and V, where each
# matrix is the discretized field of velocity values either for x or y,
# respectively.
#
##############################################################################

def give_Me_Magnitude_Velocity(U,V):
    ''' Computes velocity magnitude from two matrices of velocity values
    
    Args:
        U: x-directed velocity
        V: y-directed velocity
        
    Returns:
        uMag: Magnitude of velocity'''

    # Compute magnitude of velocity
    uMag = np.sqrt( U**2 + V**2 )

    return uMag
