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
        
    pass