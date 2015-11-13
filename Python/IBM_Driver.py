'''-------------------------------------------------------------------------

 IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
 	fluid-structure interaction models. This version of the code is based off of
	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.

 Author: Nicholas A. Battista
 Email:  nick.battista@unc.edu
 Date Created: May 27th, 2015\
 Python 3.5 port by: Christopher Strickland
 Institution: UNC-CH

 This code is capable of creating Lagrangian Structures using:
 	1. Springs
 	2. Beams (*torsional springs)
 	3. Target Points
	4. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-Tension)")

 One is able to update those Lagrangian Structure parameters, e.g., 
 spring constants, resting lengths, etc
 
 There are a number of built in Examples, mostly used for teaching purposes. 
 
 If you would like us to add a specific muscle model, 
 please let Nick (nick.battista@unc.edu) know.

----------------------------------------------------------------------------'''

import numpy as np

###############################################################################
#
# FUNCTION: Actual DRIVER of the code, where the time-stepping occurs ->
#           gets called by main2d to do the "magic" :)
#
###############################################################################

def main(struct_name, mu, rho, grid_Info, dt, T_FINAL, model_Info):

    ''' 2D IMMERSED BOUNDARY SOLVER ON RECTANGULAR DOMAIN w/ PERIODIC BOUNDARIES

    Args:
        struct_name: structure name
        mu: dynamic viscosity
        rho: density
        grid_Info: list of grid properties
        dt: time-step
        T_FINAL: final simulation time
        model_Info: list of model structure properties
        
    Returns:
        X: here
        Y: are
        U: some
        V: variables
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

    # Temporal Information
    NTime = np.floor(T_FINAL/dt)+1 # number of total time-steps,
                                # (floored, so exact number of time-steps)
    dt = T_FINAL/NTime #time-step (slightly perturbed dt, so exact number of 
                       #time-steps are used
    current_time = 0.0
    
    # GRID INFO #
    Nx = grid_Info[0]   # num of Eulerian pts. in x-direction
    Ny = grid_Info[1]   # num of Eulerian pts. in y-direction
    Lx = grid_Info[2]   # Length of Eulerian grid in x-coordinate
    Ly = grid_Info[3]   # Length of Eulerian grid in y-coordinate
    dx = grid_Info[4]   # Spatial-size in x
    dy = grid_Info[5]   # Spatial-size in y
    supp = grid_Info[6] # Delta-function support
    
    # PRINTING/PLOTTING INFO #
    pDump = grid_Info[7]        # Print (Plot) Dump interval
    pMatplotlib = grid_Info[8]  # Plot in matplotlib? (1=YES,0=NO)
    lagPlot = grid_Info[9]      # Plot LAGRANGIAN PTs ONLY in matplotlib
    velPlot = grid_Info[10]     # Plot LAGRANGIAN PTs + VELOCITY FIELD in matplotlib
    vortPlot = grid_Info[11]    # Plot LAGRANGIAN PTs + VORTICITY colormap in matplotlib
    uMagPlot = grid_Info[12]    # Plot LAGRANGIAN PTs + MAGNITUDE OF VELOCITY
                                #   colormap in matplotlib
    pressPlot = grid_Info[13]   # Plot LAGRANTIAN PTs + PRESSURE colormap in matplotlib
    
    
    # MODEL STRUCTURE DATA STORED #
    springs_Yes = model_Info[0]         # Springs: 0 (for no) or 1 (for yes) 
    update_Springs_Flag = model_Info[1] # Update_Springs: 0 (for no) or 1 (for yes)
    target_pts_Yes = model_Info[2]      # Target_Pts: 0 (for no) or 1 (for yes)
    update_Target_Pts = model_Info[3]   # Update_Target_Pts: 0 (for no) or 1 (for yes)
    beams_Yes = model_Info[4]           # Beams: 0 (for no) or 1 (for yes)
    update_Beams_Flag = model_Info[5]   # Update_Beams: 0 (for no) or 1 (for yes)
    muscles_Yes = model_Info[6]         # FV-LT Muscles: 0 (for no) or 1 (for yes)
    hill_3_muscles_Yes = model_Info[7]  # Hill 3-Element Muscle: 0 (for no) or 1 (for yes)
    arb_ext_force_Yes = model_Info[8]   # Arbitrary External Force: 0 (for no) or 1 (for yes)
    tracers_Yes = model_Info[9]         # Tracers: 0 (for no) or 1 (for yes)
    mass_Yes = model_Info[10]           # Mass Points: 0 (for no) or 1 (for yes)
    gravity_Yes = model_Info[11]        # Gravity: 0 (for no) or 1 (for yes)
    #NOTE: model_Info[12]/[13] - components of gravity vector
    porous_Yes = model_Info[14]         # Porous Media: 0 (for no) or 1 (for yes)
    concentration_Yes = model_Info[15]  # Background Concentration Gradient: 
                                        #   0 (for no) or 1 (for yes)
    
    
    
    
    #Lagrangian Structure Data
    ds = Lx/(2.*Nx)             #Lagrangian Spacing
    grid_Info[8] = ds
    
    
    # Create EULERIAN Mesh (these assume periodicity in x and y)
    x = np.arange(0,Lx,dx)
    y = np.arange(0,Ly,dy)
    # Create x-Mesh
    X = np.empty((Nx,len(x)))
    for ii in range(int(Nx)):
        X[ii,] = x
    # Create y-Mesh
    Y = np.empty((len(y),Ny))
    for ii in range(int(Ny)):
        Y[:,ii] = y
        
    # # # # # HOPEFULLY WHERE I CAN READ IN INFO!!! # # # # #


    # READ IN LAGRANGIAN POINTS #
    [Nb,xLag,yLag] = read_Vertex_Points(struct_name)
    grid_Info[7] = Nb          # num Total Number of Lagrangian Pts.
    xLag_P = xLag              # Initialize previous Lagrangian x-Values 
                               #   (for use in muscle-model)
    yLag_P = yLag              # Initialize previous Lagrangian y-Values 
                               #   (for use in muscle-model)
                            
                            
    # READ IN TRACERS (IF THERE ARE TRACERS) #
    if (tracers_Yes == 1):
       nulvar,xT,yT = read_Tracer_Points(struct_name)
       tracers = np.zeros((len(xT),4))
       tracers[0,0] = 1
       tracers[:,1] = xT
       tracers[:,2] = yT
            #tracers_info: col 1: xPt of Tracers
            #              col 2: yPt of Tracers
    else:
       tracers = 0
       
       
    # READ IN CONCENTRATION (IF THERE IS A BACKGROUND CONCENTRATION) #
    if ( concentration_Yes == 1 ):
        C,kDiffusion = read_In_Concentration_Info(struct_name)
            #C:           Initial background concentration
            #kDiffusion:  Diffusion constant for Advection-Diffusion
    else:
        C = 0 # placeholder for plotting 
        
        
    # READ IN SPRINGS (IF THERE ARE SPRINGS) #
    if ( springs_Yes == 1 ):
        springs_info = read_Spring_Points(struct_name)
            #springs_info: col 1: starting spring pt (by lag. discretization)
            #              col 2: ending spring pt. (by lag. discretization)
            #              col 3: spring stiffness
            #              col 4: spring resting lengths
    else:
        springs_info = 0  #just to pass placeholder into 
            # "please_Find_Lagrangian_Forces_On_Eulerian_grid function"
    
    
    
    
    # READ IN MUSCLES (IF THERE ARE MUSCLES) #
    if ( muscles_Yes == 1 ):
        muscles_info = read_Muscle_Points(struct_name)
            #         muscles: col 1: MASTER NODE (by lag. discretization)
            #         col 2: SLAVE NODE (by lag. discretization)
            #         col 3: length for max. muscle tension
            #         col 4: muscle constant
            #         col 5: hill parameter, a
            #         col 6: hill parameters, b
            #         col 7: force maximum!
    else:
        muscles_info = 0  #just to pass placeholder into 
            # "please_Find_Lagrangian_Forces_On_Eulerian_grid function"

    
    
    
    
    
    # READ IN MUSCLES (IF THERE ARE MUSCLES) #
    if ( hill_3_muscles_Yes == 1 ):
        muscles3_info = read_Hill_3Muscle_Points(struct_name)
            #         muscles: col 1: MASTER NODE (by lag. discretization)
            #         col 2: SLAVE NODE (by lag. discretization)
            #         col 3: length for max. muscle tension
            #         col 4: muscle constant
            #         col 5: hill parameter, a
            #         col 6: hill parameters, b
            #         col 7: force maximum!
    else:
        muscles3_info = 0  #just to pass placeholder into "please_Find_Lagrangian_Forces_On_Eulerian_grid function"
    
    
    
    
    
    # READ IN MASS POINTS (IF THERE ARE MASS PTS) #
    if ( mass_Yes == 1):
        mass_aux = read_Mass_Points(struct_name)
        #target_aux: col 0: Lag Pt. ID w/ Associated Mass Pt.
        #            col 1: "Mass-spring" stiffness parameter
        #            col 2: "MASS" value parameter
        
        # initialize mass_info
        mass_info = np.empty((len(mass_aux[:,0]),5))
        
        mass_info[:,0] = mass_aux[:,0] #Stores Lag-Pt IDs in col vector
        
        for ii in range(len(mass_info[:,0])):
            id = mass_info[ii,0]
            #here, i'm going to guess that mass pt. IDs start at 1 given prev. code
            mass_info[ii,1] = xLag[int(id)-1]  #Stores Original x-Lags as x-Mass Pt. Identities
            mass_info[ii,2] = yLag[int(id)-1]  #Stores Original y-Lags as y-Mass Pt. Identities
       
        mass_info[:,3] = mass_aux[:,1]   #Stores "mass-spring" parameter 
        mass_info[:,4] = mass_aux[:,2]   #Stores "MASS" value parameter
        
    else:
        mass_info = 0;





    # READ IN TARGET POINTS (IF THERE ARE TARGET PTS) #
    if ( target_pts_Yes == 1):
        target_aux = read_Target_Points(struct_name)
        #target_aux: col 0: Lag Pt. ID w/ Associated Target Pt.
        #            col 1: target STIFFNESSES
        
        # initialize target_info
        target_info = np.empty((len(target_aux[:,0]),4))
        
        target_info[:,0] = target_aux[:,0] #Stores Lag-Pt IDs in col vector
        for ii in range(len(target_info[:,0])):
            id = target_info[ii,0]
            #here, i'm going to guess that mass pt. IDs start at 1 given prev. code
            target_info[ii,1] = xLag[int(id)-1] #Stores Original x-Lags as 
                                                #  x-Target Pt. Identities
            target_info[ii,2] = yLag[int(id)-1] #Stores Original y-Lags as 
                                                #  y-Target Pt. Identities
       
        target_info[:,3] = target_aux[:,1] #Stores Target Stiffnesses 
    else:
        target_info = 0
    
    
    
    # READ IN POROUS MEDIA INFO (IF THERE IS POROSITY) #
    if ( porous_Yes == 1):
        porous_aux = read_Porous_Points(struct_name)
        #porous_aux: col 1: Lag Pt. ID w/ Associated Porous Pt.
        #            col 2: Porosity coefficient
        
        # initizlize porous_info
        porous_info = np.empty((len(porous_aux),4))
        
        porous_info[:,0] = porous_aux[:,0] #Stores Lag-Pt IDs in col vector
        for ii in range(porous_info[:,0]):
            id = porous_info[ii,0]
            #here, i'm going to guess that mass pt. IDs start at 1 given prev. code
            porous_info[ii,1] = xLag[int(id)-1] #Stores Original x-Lags as 
                                                #    x-Porous Pt. Identities
            porous_info[ii,2] = yLag[int(id)-1] #Stores Original y-Lags as 
                                                #    y-Porous Pt. Identities
        
        porous_info[:,3] = porous_aux[:,1] #Stores Porosity Coefficient 
    else:
        porous_info = 0



    # READ IN BEAMS (IF THERE ARE BEAMS) #
    if ( beams_Yes == 1):
        beams_info = read_Beam_Points(struct_name)
        #beams:      col 1: 1ST PT.
        #            col 2: MIDDLE PT. (where force is exerted)
        #            col 3: 3RD PT.
        #            col 4: beam stiffness
        #            col 5: curavture
    else:
        beams_info = 0


    pass
    
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

def read_In_Concentration_Info(struct_name): #untested
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
# FUNCTION: Reads in the # of springs and all MASTER NODEs, SLAVE NODEs,
#           spring STIFFNESSES, spring RESTING LENGTHS
#
###########################################################################

def read_Spring_Points(struct_name): #untested
    ''' Reads in the num of springs, master/slave nodes, spring stiffness/resting lengths
    
    Args:
        struct_name: structure name
        
    Returns:
        springs: above info stored in columns'''

    filename = struct_name+'.spring'  #Name of file to read in
    with open(filename) as f:
    #Store elements on .spring file into a matrix starting w/ 2nd row of read in data.
        springs = np.loadtxt(f,skiprows=1,usecols=(0,1,2,3))

    #springs: col 1: starting spring pt (by lag. discretization)
    #         col 2: ending spring pt. (by lag. discretization)
    #         col 3: spring stiffness
    #         col 4: spring resting lengths
    
    return springs



    
###########################################################################
#
# FUNCTION: Reads in the # of muscles and all MASTER NODEs, SLAVE NODEs,
#           length for max. muscle tension, muscle constant, hill
#           parameters (a and b), and Force-Max
#
###########################################################################

def read_Muscle_Points(struct_name): #untested
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

def read_Mass_Points(struct_name): #untested
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

def read_Target_Points(struct_name): #untested
    ''' Reads in the num of target pts, target-pt-nodes, and target-stiffness
    
    Args:
        struct_name: structure name
        
    Returns:
        targets: array of target info'''
    filename = struct_name+'.target'  #Name of file to read in
    with open(filename) as f:
        targets = np.txtload(f,skiprows=1,usecols=(0,1))

    #targets: col 1: Lag Pt. ID w/ Associated Target Pt.
    #         col 2: target STIFFNESSES
    
    return targets
    

    
###########################################################################
#
# FUNCTION: Reads in the # of POROUS PTS, POROUS-PT-NODEs, and their
#           POROUSITY-COEFFICIENTS
#
###########################################################################

def read_Porous_Points(struct_name): #untested
    ''' Reads in the num of porous pts, pt-nodes, and porousity coefficients
    
    Args:
        struct_name: structure name
        
    Returns:
        porosity: array of porosity info'''

    filename = struct_name+'.porous'  #Name of file to read in
    with open(filename) as f:
        porosity = np.loadtxt(f,skiprows=1,usecols=(0,1))

    #porous:  col 1: Lag Pt. ID w/ Associated Porous Pt.
    #         col 2: Porosity coefficient
    
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