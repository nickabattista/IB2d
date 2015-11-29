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

import pdb
import numpy as np
from math import sqrt
import os
from Supp import *
from please_Find_Lagrangian_Forces_On_Eulerian_grid import\
    please_Find_Lagrangian_Forces_On_Eulerian_grid

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

    
    # Temporal Information
    NTime = np.floor(T_FINAL/dt)+1 # number of total time-steps,
                                # (floored, so exact number of time-steps)
    dt = T_FINAL/NTime #time-step (slightly perturbed dt, so exact number of 
                       #time-steps are used
    current_time = 0.0
    
    # GRID INFO #
    Nx = int(grid_Info[0])   # num of Eulerian pts. in x-direction
    Ny = int(grid_Info[1])   # num of Eulerian pts. in y-direction
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
    X = np.empty((Nx,x.size))
    for ii in range(Nx):
        X[ii,] = x
    # Create y-Mesh
    Y = np.empty((y.size,Ny))
    for ii in range(Ny):
        Y[:,ii] = y
        
    # # # # # HOPEFULLY WHERE I CAN READ IN INFO!!! # # # # #

    # READ IN LAGRANGIAN POINTS #
    Nb,xLag,yLag = read_Vertex_Points(struct_name)
    grid_Info[7] = Nb          # num Total Number of Lagrangian Pts.
    xLag_P = xLag              # Initialize previous Lagrangian x-Values 
                               #   (for use in muscle-model)
    yLag_P = yLag              # Initialize previous Lagrangian y-Values 
                               #   (for use in muscle-model)
                            
                            
    # READ IN TRACERS (IF THERE ARE TRACERS) #
    if (tracers_Yes == 1):
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
        springs_info = np.zeros((1,1))  #just to pass placeholder into 
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
        muscles_info = np.zeros((1,1))  #just to pass placeholder into 
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
        muscles3_info = np.zeros((1,1))  #just to pass placeholder into "please_Find_Lagrangian_Forces_On_Eulerian_grid function"
    
    
    
    
    
    # READ IN MASS POINTS (IF THERE ARE MASS PTS) #
    if ( mass_Yes == 1):
        mass_aux = read_Mass_Points(struct_name)
        #target_aux: col 0: Lag Pt. ID w/ Associated Mass Pt.
        #            col 1: "Mass-spring" stiffness parameter
        #            col 2: "MASS" value parameter
        
        # initialize mass_info
        mass_info = np.empty((mass_aux.shape[0],5))
        
        mass_info[:,0] = mass_aux[:,0] #Stores Lag-Pt IDs in col vector
        
        for ii in range(mass_info.shape[0]):
            id = int(mass_info[ii,0])
            #here, i'm going to guess that mass pt. IDs start at 1 given prev. code
            mass_info[ii,1] = xLag[id-1]  #Stores Original x-Lags as x-Mass Pt. Identities
            mass_info[ii,2] = yLag[id-1]  #Stores Original y-Lags as y-Mass Pt. Identities
       
        mass_info[:,3] = mass_aux[:,1]   #Stores "mass-spring" parameter 
        mass_info[:,4] = mass_aux[:,2]   #Stores "MASS" value parameter
        
    else:
        mass_info = np.zeros((1,1))





    # READ IN TARGET POINTS (IF THERE ARE TARGET PTS) #
    if ( target_pts_Yes == 1):
        target_aux = read_Target_Points(struct_name)
        #target_aux: col 0: Lag Pt. ID w/ Associated Target Pt.
        #            col 1: target STIFFNESSES
        
        # initialize target_info
        target_info = np.empty((target_aux.shape[0],4))
        
        target_info[:,0] = target_aux[:,0] #Stores Lag-Pt IDs in col vector
        for ii in range(target_info.shape[0]):
            id = int(target_info[ii,0])
            #here, i'm going to guess that mass pt. IDs start at 1 given prev. code
            target_info[ii,1] = xLag[id-1] #Stores Original x-Lags as 
                                                #  x-Target Pt. Identities
            target_info[ii,2] = yLag[id-1] #Stores Original y-Lags as 
                                                #  y-Target Pt. Identities
       
        target_info[:,3] = target_aux[:,1] #Stores Target Stiffnesses 
    else:
        target_info = np.zeros((1,1))
    
    
    
    # READ IN POROUS MEDIA INFO (IF THERE IS POROSITY) #
    if ( porous_Yes == 1):
        porous_aux = read_Porous_Points(struct_name)
        #porous_aux: col 1: Lag Pt. ID w/ Associated Porous Pt.
        #            col 2: Porosity coefficient
        
        # initizlize porous_info
        porous_info = np.empty((porous_aux.size,4))
        
        porous_info[:,0] = porous_aux[:,0] #Stores Lag-Pt IDs in col vector
        for ii in range(porous_info.shape[0]):
            id = int(porous_info[ii,0])
            #here, i'm going to guess that mass pt. IDs start at 1 given prev. code
            porous_info[ii,1] = xLag[id-1] #Stores Original x-Lags as 
                                                #    x-Porous Pt. Identities
            porous_info[ii,2] = yLag[id-1] #Stores Original y-Lags as 
                                                #    y-Porous Pt. Identities
        
        porous_info[:,3] = porous_aux[:,1] #Stores Porosity Coefficient 
    else:
        porous_info = np.zeros((1,1))



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
        
    
    # CONSTRUCT GRAVITY INFORMATION (IF THERE IS GRAVITY) #
    if gravity_Yes == 1:
        #gravity_Vec[0] = model_Info[11]     # x-Component of Gravity Vector
        #gravity_Vec[1] = model_Info[12]     # y-Component of Gravity Vector
        xG = model_Info[12]
        yG = model_Info[13]
        normG = sqrt( xG**2 + yG**2 )
        gravity_Info = [gravity_Yes, xG/normG, yG/normG]
        #   col 1: flag if considering gravity
        #   col 2: x-component of gravity vector (normalized)
        #   col 3: y-component of gravity vector (normalized)
        
        del xG, yG, normG
        
    else:
        gravity_Info = np.zeros((1,1))

    
    # Initialize the initial velocities to zero.
    U = np.zeros((Ny,Nx))                           # x-Eulerian grid velocity
    V = U                                           # y-Eulerian grid velocity
    mVelocity = np.zeros((mass_info.shape[0],2))  # mass-Pt velocity 

    if arb_ext_force_Yes == 1:
        firstExtForce = 1                           # initialize external forcing
        indsExtForce = 0                            # initialize for external forcing computation
    
    # ACTUAL TIME-STEPPING IBM SCHEME! 
    #(flags for storing structure connects for printing and printing to .vtk)
    cter = 0; ctsave = 0; firstPrint = 1; loc = 1; diffy = 1
    
    # CREATE VIZ_IB2D FOLDER and VISIT FILES
    try:
        os.mkdir('viz_IB2d')
    except FileExistsError:
        #File already exists
        pass
    #I'm going to expect that vizID is a file object with write permission...?
    vizID = 1 #JUST INITIALIZE BC dumps.visit isn't working correctly...yet
    os.chdir('viz_IB2d')
    #vizID = open('dumps.visit','w')
    #vizID.write('!NBLOCKS 6\n')
    #os.chdir('..')
    
    #Initialize Vorticity, uMagnitude, and Pressure for initial colormap
    #Print initializations to .vtk
    vort = np.zeros((Ny,Nx)); uMag = np.array(vort); p = np.array(vort)
    lagPts = np.zeros((xLag.size,3))
    lagPts[:,0] = xLag; lagPts[:,1] = yLag
    connectsMat,spacing = give_Me_Lag_Pt_Connects(ds,xLag,yLag,Nx)
    print_vtk_files(ctsave,vizID,vort,uMag,p,U,V,Lx,Ly,Nx,Ny,lagPts,\
    connectsMat,tracers,concentration_Yes,C)
    print('Current Time(s): {0}\n'.format(current_time))
    ctsave += 1
    
    #
    #
    # * * * * * * * * * * BEGIN TIME-STEPPING! * * * * * * * * * * *
    #
    #
    
    #This line commented for debugging
    #while current_time < T_FINAL:
        
    #
    #******Step 1: Update Position of Boundary of membrane at half time-step ******
    #                 (Variables end with h if it is a half-step)
    #
    xLag_h,yLag_h = please_Move_Lagrangian_Point_Positions(U, V, xLag, yLag,\
        xLag, yLag, x, y, dt/2, grid_Info, 0)
        
    if mass_Yes == 1:
        mass_info, massLagsOld = please_Move_Massive_Boundary(dt/2,\
        mass_info,mVelocity)
       
    if ( ( update_Springs_Flag == 1 ) and ( springs_Yes == 1 ) ):
        #This function is application specific, located with main2d
        springs_info = update_Springs(dt,current_time,xLag,yLag,springs_info)
        
    if ( ( update_Target_Pts == 1 ) and ( target_pts_Yes == 1) ):
        #This function is application specific, located with main2d
        target_info = update_Target_Point_Positions(dt,current_time,target_info)
        
    if ( ( update_Beams_Flag == 1 ) and ( beams_Yes == 1) ):
        #This function is application specific, located with main2d
        beams_info = update_Beams(dt,current_time,beams_info)
        
    #
    #*******STEP 2: Calculate Force coming from membrane at half time-step ********
    #
    Fxh, Fyh, F_Mass_Bnd, F_Lag = please_Find_Lagrangian_Forces_On_Eulerian_grid(\
    dt, current_time, xLag_h, yLag_h, xLag_P, yLag_P, x, y, grid_Info, model_Info,\
    springs_info, target_info, beams_info, muscles_info, muscles3_info, mass_info)

    
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
    
    
##############################################################################
#
# FUNCTION: give me Connects Vector for printing Lagrangian .vtk info!
#
##############################################################################

def give_Me_Lag_Pt_Connects(ds,xLag,yLag,Nx):
    ''' Give me Connects Vector for printing Lagrangian .vtk info!

    Args:
        ds:
        xLag:
        yLag:
        Nx:
        
    Returns:
        connectsMat:
        space:'''

    N = xLag.size

    if Nx <= 32:
        space = 5*ds
    elif Nx <= 64:
       space = 5*ds
    elif Nx <=128:
       space = 5*ds
    elif Nx <=256:
        space = 10*ds
    elif Nx <= 512:
        space = 20*ds
    else:
        space = 40*ds
        

    #need to instantiate connectsMat or something...
    connectsMat0 = []; connectsMat1 = []
    for ii in range(N): #for i=1:N
        if ii<N-1:
            x1=xLag[ii]; x2=xLag[ii+1]
            y1=yLag[ii]; y2=yLag[ii+1]
            dist = sqrt( (x1-x2)**2 + (y1-y2)**2 )
            if dist < space:
                #The note here refers to Cpp notation (and .vtk counting)...
                #   I'm guessing that means counting starts at 0 as in Python
                connectsMat0.append(ii)   #For Cpp notation (and .vtk counting)
                connectsMat1.append(ii+1) #For Cpp notation (and .vtk counting)
        elif ii==N-1:
            x1=xLag[ii]; x2=xLag[0]
            y1=yLag[ii]; y2=yLag[0]
            dist = sqrt( (x1-x2)**2 + (y1-y2)**2 )
            if dist < space:
                connectsMat0.append(N-1) #For Cpp notation (and .vtk counting)
                connectsMat1.append(0)   #For Cpp notation (and .vtk counting)
    connectsMat = np.array([connectsMat0,connectsMat1]).T
    
    return (connectsMat,space)


##############################################################################
#
# FUNCTION: gives appropriate string number for filename in printing the
# .vtk files.
#
##############################################################################

def print_vtk_files(ctsave,vizID,vort,uMag,p,U,V,Lx,Ly,Nx,Ny,lagPts,\
    connectsMat,tracers,concentration_Yes,C):
    ''' Gives appropriate string number for filename in printing the .vtk files'''

    #Give spacing for grid
    dx = Lx/Nx 
    dy = Ly/Ny


    #Go into viz_IB2d directory. This was throwing an error because we're already there!
    if os.path.split(os.getcwd())[1] != 'viz_IB2d':
        os.chdir('viz_IB2d')

    #Find string number for storing files
    strNUM = give_String_Number_For_VTK(ctsave)
    vortfName = 'Omega.'+strNUM+'.vtk'
    uMagfName = 'uMag.'+strNUM+'.vtk'
    pfName = 'P.'+strNUM+'.vtk'
    uXName = 'uX.'+strNUM+'.vtk'
    uYName = 'uY.'+strNUM+'.vtk'
    velocityName = 'u.'+strNUM+'.vtk'
    lagPtsName = 'lagsPts.'+strNUM+'.vtk'
    lagPtsConName = 'lagPtsConnect.'+strNUM+'.vtk'

    #Print Lagrangian Pts to .vtk format
    savevtk_points(lagPts, lagPtsName, 'lagPts')

    #Print Lagrangian Pts w/ CONNECTIONS to .vtk format
    savevtk_points_connects(lagPts, lagPtsConName, 'lagPtsConnected',connectsMat)

    #Print Tracer Pts (*if tracers*)
    if tracers[0,0] == 1:
        tracersPtsName = 'tracer.'+strNUM+'.vtk'
        #tMatrix = tracers[:,1:4]
        savevtk_points(tracers[:,1:4],tracersPtsName, 'tracers') 
            
    #Print another cycle to .visit file
    #vizID.write(vortfName+'\n')
    #vizID.write(uMagfName+'\n')
    #vizID.write(pfName+'\n')
    #vizID.write(uXName+'\n')
    #vizID.write(uYName+'\n')
    #vizID.write(velocityName+'\n')


    #Print SCALAR DATA (i.e., colormap data) to .vtk file
    savevtk_scalar(vort, vortfName, 'Omega',dx,dy)
    savevtk_scalar(uMag, uMagfName, 'uMag',dx,dy)
    savevtk_scalar(p, pfName, 'P',dx,dy)
    savevtk_scalar(U, uXName, 'uX',dx,dy)
    savevtk_scalar(V, uYName, 'uY',dx,dy)

    if concentration_Yes == 1:
        confName = 'concentration.'+strNUM+'.vtk'
        savevtk_scalar(C, confName, 'Concentration',dx,dy)

    #Print VECTOR DATA (i.e., velocity data) to .vtk file
    savevtk_vector(U, V, velocityName, 'u',dx,dy)

    #Get out of viz_IB2d folder
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

    #TRY PRINTING THEM AS UNSTRUCTURED_GRID
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
            file.write('{0} {1:d} {2:d}\n'.format(2,connectsMat[s,0],connectsMat[s,1]))
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


    #TRY PRINTING THEM AS UNSTRUCTURED_GRID
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
        file.write('CELLS {0} {1}\n'.format(N,2*N))
        for s in range(N):
            file.write('{0} {1}\n'.format(1,s))
        file.write('\n')
        #
        file.write('CELL_TYPES {0}\n'.format(N)) # N = # of "Cells"
        for ii in range(N):
           file.write('1 ')
        file.write('\n')



    #TRY PRINTING THEM AS POLYGONAL DATA
    # with open(filename,'w') as file:
        # file.write('# vtk DataFile Version 2.0\n')
        # file.write(vectorName+'\n')
        # file.write('ASCII\n')
        # file.write('DATASET STRUCTURED_GRID\n')
        # file.write('DIMENSIONS 64 1 1\n')
        # file.write('POINTS {0} float\n', N)
        # for ii in range(N):
            # file.write('{0:.15e} {1:.15e} {2:.15e}\n'.format(X[ii,0],X[ii,1],X[ii,2]))
        # file.write('1.1 1.1 0\n')
        # file.write('CELL_DATA 1\n')
        # file.write('POINT_DATA {0} \n',N)
        # file.write('FIELD FieldData 1\n')
        # file.write('nodal 1 {0} float\n'.format(N)
        # file.write('0 1 1.1 2\n')
        # file.write('SCALARS nodal float\n')
        # file.write('SCALARS '+vectorName+' float 1 \n')
        # file.write('LOOKUP_TABLE default\n')


    # TRY PRINTING THEM AS POINTS
    # with open(filename,'w') as file:
        # file.write('# vtk DataFile Version 2.0\n')
        # file.write('Cube example\n')
        # file.write('ASCII\n')
        # file.write('DATASET UNSTRUCTURED_GRID\n')
        # file.write('POINTS {0} float\n'.format(N))
        # for ii in range(N):
            # file.write('{0:.15e} {1:.15e} {2:.15e}\n'.format(X[ii,0],X[ii,1],X[ii,2]))
        # file.write('POINT_DATA {0} \n'.format(N)
        # file.write('SCALARS '+vectorName+' float 1 \n')
        # file.write('LOOKUP_TABLE default\n')



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
        vectorName:
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
    #   So, specifically, nz is now gone. I will keep the output the same,
    #   however, for compatibility. So 1 will be pritned in the Z column.
    
    #I'm changing the fprintf that was here to an error. If you want, you can
    #   catch it in the following function.
    assert (X.shape == Y.shape), 'Error: velocity arrays of unequal size'
    nx, ny = X.shape
    #Python 3.5 automatically opens in text mode unless otherwise specified
    with open(filename,'w') as fid:
        fid.write('# vtk DataFile Version 2.0\n')
        fid.write('Comment goes here\n')
        fid.write('ASCII\n')
        fid.write('\n')
        fid.write('DATASET STRUCTURED_POINTS\n')
        # 1 below was nz
        fid.write('DIMENSIONS    {0}   {1}   {2}\n'.format(nx, ny, 1))
        fid.write('\n')
        fid.write('ORIGIN    0.000   0.000   0.000\n')
        #fid.write('SPACING   1.000   1.000   1.000\n') #if want [1,32]x[1,32] rather than [0,Lx]x[0,Ly]
        fid.write('SPACING   '+str(dx)+str(dy)+'   1.000\n')
        fid.write('\n')
        fid.write('POINT_DATA   {0}\n'.format(nx*ny))
        fid.write('VECTORS '+vectorName+' double\n')
        fid.write('\n')
        for b in range(ny):
            for c in range(nx):
                fid.write('{0} '.format(X[c,b]))
                fid.write('{0} '.format(Y[c,b]))
                fid.write('1 ')
            fid.write('\n')


##############################################################################
#
# FUNCTION: prints scalar matrix to vtk formated file
#
##############################################################################

def savevtk_scalar(array, filename, colorMap,dx,dy):
    ''' Prints scalar matrix to vtk formatted file.
    
    Args:
        array: 2-D ndarray
        filename: file name
        colorMap:
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
    #   So, specifically, nz is now gone. I will keep the output the same,
    #   however, for compatibility. So 1 will be pritned in the Z column.
    nx,ny = array.shape
    #Python 3.5 automatically opens in text mode unless otherwise specified
    with open(filename,'w') as fid:
        fid.write('# vtk DataFile Version 2.0\n')
        fid.write('Comment goes here\n')
        fid.write('ASCII\n')
        fid.write('\n')
        fid.write('DATASET STRUCTURED_POINTS\n')
        # 1 below was nz
        fid.write('DIMENSIONS    {0}   {1}   {2}\n'.format(nx, ny, 1))
        fid.write('\n')
        fid.write('ORIGIN    0.000   0.000   0.000\n')
        #fid.write('SPACING   1.000   1.000   1.000\n') #if want [1,32]x[1,32] rather than [0,Lx]x[0,Ly]
        fid.write('SPACING   '+str(dx)+str(dy)+'   1.000\n')
        fid.write('\n')
        # The 1 below was nz
        fid.write('POINT_DATA   {0}\n'.format(nx*ny*1))
        fid.write('SCALARS '+colorMap+' double\n')
        fid.write('LOOKUP_TABLE default\n')
        fid.write('\n')
        for b in range(ny):
            for c in range(nx):
                fid.write('{0} '.format(array[c,b]))
            fid.write('\n')