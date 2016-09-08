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

from math import sqrt
import numpy as np
from Supp import give_1D_NonZero_Delta_Indices
from Supp import give_Eulerian_Lagrangian_Distance, give_Delta_Kernel
 
################################################################################
#
# FUNCTION: Computes the components of the force term in Navier-Stokes from
#           deformations of the boundary of the immersed boundary
#
################################################################################


def please_Find_Lagrangian_Forces_On_Eulerian_grid(dt, current_time, xLag, yLag,\
    xLag_P,yLag_P, x, y, grid_Info, model_Info, springs, targets, beams,\
    muscles, muscles3, masses, d_Springs):
    ''' Compute components of force term in Navier-Stokes from boundary deformations

        Args:
            dt:
            current_time: Current time of simulation (in seconds)
            xLag:         x positions of Lagrangian structure
            yLag:         y positions of Lagrangian structure
            xLag_P:
            yLag_P:
            x:            x positions on Eulerian grid
            y:            y positions on Eulerian grid
            grid_Info:    holds lots of geometric pieces about grid / simulations
            model_Info:   Stores if springs, if update_springs, if target_pts, if 
                update_target_pts (as 0 (no) or 1 (yes) )
            springs:      Stores Master Node, Slave Node, Spring Stiffness, 
                Restling-Lengths, all in column vecs
            targets:      Stores target point index, correponding xLag, yLag and 
                target point stiffness
            beams:        Stores 1st Node, 2nd (MIDDLE-MAIN) Node, 3rd Nodes, 
                Beam Stiffnesses, and Beam curvatures
            muscles:
            muscles3:
            masses:       Stores mass point index, correponding xLag, yLag, 
                "spring" stiffness, and mass value parameter
                
        Returns:
            Fx:
            Fy:
            F_Mass:
            F_Lag:
            
    The components of the force are given by
    F(x,y) = int{ f(s) * delta(x - xLag(s)) * delta(y - yLag(s)) * ds }
     
    where s parameteriizes the Lagrangian structure.

    Force density is computed using a SIMPLE LINEAR SPRING model w/ resting length L.  
    This leads to a force density of the form,
                   f = ( \LagPts_s * ( 1 - L / abs(\LagPts_s) )  )/ds^2'''

    # Grid Info - list#
    Nx =    grid_Info['Nx'] # # of Eulerian pts. in x-direction
    Ny =    grid_Info['Ny'] # # of Eulerian pts. in y-direction
    Lx =    grid_Info['Lx'] # Length of Eulerian grid in x-coordinate
    Ly =    grid_Info['Ly'] # Length of Eulerian grid in y-coordinate
    dx =    grid_Info['dx'] # Spatial-size in x
    dy =    grid_Info['dy'] # Spatial-size in y
    supp =  grid_Info['supp'] # Delta-function support
    Nb =    grid_Info['Nb'] # # of Lagrangian pts. 
    ds =    grid_Info['ds'] # Lagrangian spacing


    # Model Potential Forces #
    springs_Yes = model_Info['springs']               # Springs: 0 (for no) or 1 (for yes) 
    target_pts_Yes = model_Info['target_pts']         # Target_Pts: 0 (for no) or 1 (for yes)
    beams_Yes = model_Info['beams']                   # Beams: 0 (for no) or 1 (for yes)
    muscle_LT_FV_Yes = model_Info['muscles']          # Length-Tension/Force-Velocity Muscle: 
        # 0 (for no) or 1 (for yes) (Length/Tension - Hill Model)
    muscle_3_Hill_Yes = model_Info['hill_3_muscles']  # 3-Element Hill Model: 
        # 0 (for no) or 1 (for yes) (3 Element Hill + Length-Tension/Force-Velocity)
    mass_Yes = model_Info['mass']                     # Mass Pts: 0 (for no) or 1 (for yes)
    d_Springs_Yes = model_Info['damped_springs']      # Damped Springs: 0 (for no) or 1 (for yes)

    #
    # Compute MUSCLE LENGTH-TENSION/FORCE-VELOCITY #
    #(if using combined length/tension-Hill model) #
    #
    
    
    if ( muscle_LT_FV_Yes == 1):
        fx_muscles, fy_muscles = give_Muscle_Force_Densities(Nb,xLag,yLag,\
            xLag_P,yLag_P,muscles,current_time,dt)
    else:
        fx_muscles = np.zeros(xLag.size)
        fy_muscles = np.zeros(xLag.size)

    

    #
    # Compute 3-ELEMENT HILL MUSCLE MODEL FORCE DENSITIES (if using combined 3-HILL + LT/FV) #
    #
    if ( muscle_3_Hill_Yes == 1):
        fx_muscles3, fy_muscles3 = give_3_Element_Muscle_Force_Densities(Nb,\
            xLag,yLag,xLag_P,yLag_P,muscles3,current_time,dt)
    else:
        fx_muscles3 = np.zeros(xLag.size)
        fy_muscles3 = np.zeros(xLag.size)



    # Compute SPRING FORCE DENSITIES (if there are springs!)
    if ( springs_Yes == 1 ):
        
        # Compute "Connections Matrix" for what springs are attached to whom     #
        # (and "Connections Stiffness Matrix" and "Connections Restling Lengths" #
        #connects = give_Me_Spring_Connections_Matrix(Nb,Nsprings,sp_1,sp_2,K_Vec,L_Vec); 
        
        # Compute distances between Lag-Pts w/ Springs for Spring-Tension Calc.
        #dLag_x = give_Spring_Lagrangian_Distance(xLag, Lx, springs);
        #dLag_y = give_Spring_Lagrangian_Distance(yLag, Ly, springs);

        # Compute the Lagrangian SPRING tensions!
        #[Tx Ty] = give_Me_Spring_Lagrangian_Tension(Nb,dLag_x,dLag_y,springs);

        # Compute the Lagrangian SPRING force densities!
        fx_springs, fy_springs = give_Me_Spring_Lagrangian_Force_Densities(ds,Nb,\
            xLag,yLag,springs)
        
    else:
        fx_springs = np.zeros(Nb) #No x-forces coming from springs
        fy_springs = np.zeros(Nb) #No y-forces coming from springs

        

    # Compute MASS PT FORCE DENSITIES (if there are mass points!)
    if ( mass_Yes == 1):
        # Compute the Lagrangian MASS PT force densities!
        fx_mass, fy_mass, F_Mass = give_Me_Mass_Lagrangian_Force_Densities(ds,\
            xLag,yLag,masses)
    else:
        fx_mass = np.zeros(Nb) #No x-forces coming from mass points
        fy_mass = np.zeros(Nb) #No y-forces coming from mass points
        F_Mass = 0             #Dummy to pass along  



    # Compute TARGET FORCE DENSITIES (if there are target points!)
    if ( target_pts_Yes == 1):
        # Compute the Lagrangian TARGET force densities!
        fx_target, fy_target = give_Me_Target_Lagrangian_Force_Densities(ds,\
            xLag,yLag,targets)
        
    else:
        fx_target = np.zeros(Nb) #No x-forces coming from target points
        fy_target = np.zeros(Nb) #No y-forces coming from target points




    # Compute BEAM FORCE DENSITIES (if there are beams!)
    if ( beams_Yes == 1 ):

        # Compute the Lagrangian SPRING force densities!
        fx_beams, fy_beams = give_Me_Beam_Lagrangian_Force_Densities(ds,Nb,\
            xLag,yLag,beams)
        
    else:
        fx_beams = np.zeros(Nb) #No x-forces coming from beams
        fy_beams = np.zeros(Nb) #No y-forces coming from beams


    # Compute SPRING FORCE DENSITIES (if there are springs!)
    if ( d_Springs_Yes == 1 ):
     
        fx_dSprings, fy_dSprings = give_Me_Damped_Springs_Lagrangian_Force_Densities(ds,Nb,\
            xLag,yLag,d_Springs,xLag_P,yLag_P,dt)
        
    else:
        fx_dSprings = np.zeros(Nb) #No x-forces coming from damped springs
        fy_dSprings = np.zeros(Nb) #No y-forces coming from damped springs




    # SUM TOTAL FORCE DENSITY! #
    fx = fx_springs + fx_target + fx_beams + fx_muscles + fx_muscles3 + fx_mass + fx_dSprings
    fy = fy_springs + fy_target + fy_beams + fy_muscles + fy_muscles3 + fy_mass + fy_dSprings


    # SAVE LAGRANGIAN FORCES
    F_Lag = np.empty((fx.size,2))
    F_Lag[:,0] = fx
    F_Lag[:,1] = fy
        

    # Give me delta-function approximations!
    delta_X, delta_Y = give_Me_Delta_Function_Approximations_For_Force_Calc(x,y,\
    grid_Info,xLag,yLag)


    # Transform the force density vectors into diagonal matrices
    fxds = np.diag(fx*ds)
    fyds = np.diag(fy*ds)


    # Find Eulerian forces on grids by approximating the line integral, 
    #       F(x,y) = int{ f(s) delta(x - xLag(s)) delta(y - yLag(s)) ds }
    Fx = delta_Y @ fxds @ delta_X
    Fy = delta_Y @ fyds @ delta_X

    return (Fx, Fy, F_Mass, F_Lag)
    
    
    
################################################################################
#
# FUNCTION computes the Lagrangian MUSCLE Force Densities for 
#          LENGTH-TENSION/FORCE-VELOCITY MUSCLE MODEL.
#
################################################################################

def give_Muscle_Force_Densities(Nb,xLag,yLag,xLag_P,yLag_P,muscles,current_time,dt):
    ''' Compute Lagrangian muscle force densities
    
    Args:
        Nb:
        xLag:
        yLag:
        xLag_P:
        yLag_P:
        muscles: ndarray
        current_time:
        dt:
        
    Returns:
        fx:
        fy:'''
    
    # These just provide views (by reference) into slices of muscles.
    #   They are not assigned to, so aliasing is faster.
    
    Nmuscles = muscles.shape[0]      # # of Muscles
    id_Master = muscles[:,0].astype('int')     # MASTER NODE Muscle Connections
    id_Slave = muscles[:,1].astype('int')      # SLAVE NODE Muscle Connections
    LFO_Vec = muscles[:,2] # Stores length for max. muscle tension
    SK_Vec = muscles[:,3]  # Stores muscle constant
    a_Vec = muscles[:,4]   # Stores Hill Parameter, a
    b_Vec = muscles[:,5]   # Stores Hill Parameter, b
    FMAX_Vec = muscles[:,6]# Stores Force-Maximum for Muscle

    fx = np.zeros(Nb)                 # Initialize storage for x-forces
    fy = np.zeros(Nb)                 # Initialize storage for y-forces
    
    xPt = xLag[id_Master]        # x-Pt of interest at the moment to drive 
                                 #  muscle contraction
    
    dx = xLag[id_Slave] - xLag[id_Master]  # x-Distance btwn slave and master node
    dy = yLag[id_Slave] - yLag[id_Master]  # y-Distance btwn slave and master node
    LF = np.sqrt( dx**2 + dy**2 )          # Euclidean DISTANCE between master and slave node
    
    
    dx_P = xLag_P[id_Slave] - xLag_P[id_Master]  # x-Distance btwn slave and master node
    dy_P = yLag_P[id_Slave] - yLag_P[id_Master]  # y-Distance btwn slave and master node
    LF_P = np.sqrt( dx_P**2 + dy_P**2 )          # Euclidean DISTANCE between master and slave node
    
    v =  np.abs(LF-LF_P)/dt        # How fast the muscle is contracting/expanding
    
    # Find actual muscle activation magnitude
    # This function is user defined and included with main2d.
    #   In this Python code, it takes in vectors instead of scalars.
    #   Universal numpy functions could probably deal with this if you don't 
    #   want to port the vectorization.
    #   BE CAREFUL!! xLag is MUTABLE. Pass it by value instead of reference.
    from give_Muscle_Activation import give_Muscle_Activation
    Fm = give_Muscle_Activation(v,LF,LFO_Vec,SK_Vec,a_Vec,b_Vec,FMAX_Vec,\
                                current_time,xPt,np.array(xLag))
    
    mF_x = Fm*(dx/LF)           # cos(theta) = dx / LF;
    mF_y = Fm*(dy/LF)           # sin(theta) = dy / LF;

    for ii in range(Nmuscles):
        
        fx[id_Master[ii]] += mF_x[ii]  # Sum total forces for node,
                        # i in x-direction (this is MASTER node for this spring)
        fy[id_Master[ii]] += mF_y[ii]  # Sum total forces for node,
                        # i in y-direction (this is MASTER node for this spring)
        
        fx[id_Slave[ii]] -= mF_x[ii]    # Sum total forces for node,
                        # i in x-direction (this is SLAVE node for this spring)
        fy[id_Slave[ii]] -= mF_y[ii]    # Sum total forces for node,
                        # i in y-direction (this is SLAVE node for this spring)
                        
    return (fx,fy)
        
        

################################################################################
#
# FUNCTION: computes the Lagrangian MUSCLE Force Densities for 3-ELEMENT
#           HILL MODEL! ( coupled w/ force-velocity/length-tension model for
#           contractile element!)
#
################################################################################

def give_3_Element_Muscle_Force_Densities(Nb,xLag,yLag,xLag_P,yLag_P,muscles3,\
    current_time,dt):
    ''' Computes the Lagrangian MUSCLE force densities for 3-element hill model
    
    Args:
        Nb:
        xLag:
        yLag:
        xLag_P:
        yLag_P:
        muscles3:
        current_time:
        dt:
        
    Returns:
        fx:
        fy:'''

    # These just provide views (by reference) into slices of muscles3.
    #   They are not assigned to, so aliasing is faster.
        
    Nmuscles = muscles3.shape[0]       # # of Muscles
    id_Master = muscles3[:,0].astype('int')      # MASTER NODE Muscle Connections
    id_Slave = muscles3[:,1].astype('int')       # SLAVE NODE Muscle Connections
    LFO_Vec = muscles3[:,2]  # Stores length for max. muscle tension
    SK_Vec = muscles3[:,3]   # Stores muscle constant
    a_Vec = muscles3[:,4]    # Stores Hill Parameter, a
    b_Vec = muscles3[:,5]    # Stores Hill Parameter, b
    FMAX_Vec = muscles3[:,6] # Stores Force-Maximum for Muscle

    fx = np.zeros(Nb)                 # Initialize storage for x-forces
    fy = np.zeros(Nb)                 # Initialize storage for y-forces
    
    xPt = xLag[id_Master]     # x-Pt of interest at the moment to drive muscle contraction
    
    dx = xLag[id_Slave] - xLag[id_Master] # x-Distance btwn slave and master node
    dy = yLag[id_Slave] - yLag[id_Master] # y-Distance btwn slave and master node
    LF = np.sqrt( dx**2 + dy**2 )         # Euclidean DISTANCE between
    
    dx_P = xLag_P[id_Slave] - xLag_P[id_Master] # x-Distance btwn slave and master node
    dy_P = yLag_P[id_Slave] - yLag_P[id_Master] # y-Distance btwn slave and master node
    LF_P = np.sqrt( dx_P**2 + dy_P**2 )         # Euclidean DISTANCE between 
                                                #   master and slave node
                                                
    v =  np.abs(LF-LF_P)/dt         # How fast the muscle is contracting/expanding
    
    # Find actual muscle activation magnitude
    # This function is user defined and included with main2d.
    #   In this Python code, it takes in vectors instead of scalars.
    #   Universal numpy functions could probably deal with this if you don't 
    #   want to port the vectorization.
    #   BE CAREFUL!! xLag is MUTABLE. Pass it by value instead of reference.
    Fm = give_3_Element_Muscle_Activation(v,LF,LFO_Vec,SK_Vec,a_Vec,b_Vec,\
                                    FMAX_Vec,current_time,xPt,np.array(xLag))
        
    mF_x = Fm*(dx/LF)           # cos(theta) = dx / LF;
    mF_y = Fm*(dy/LF)           # sin(theta) = dy / LF;
    
    for ii in range(Nmuscles):
        
        fx[id_Master[ii]] += mF_x[ii]  # Sum total forces for node,
                        # i in x-direction (this is MASTER node for this spring)
        fy[id_Master[ii]] += mF_y[ii]  # Sum total forces for node, 
                        # i in y-direction (this is MASTER node for this spring)
        
        fx[id_Slave[ii]] -= mF_x[ii]    # Sum total forces for node, 
                        # i in x-direction (this is SLAVE node for this spring)
        fy[id_Slave[ii]] -= mF_y[ii]    # Sum total forces for node, 
                        # i in y-direction (this is SLAVE node for this spring)
        
    return (fx,fy)
    
    
    
###########################################################################
#
# FUNCTION computes the Lagrangian SPRING Force Densities .
#
###########################################################################

def give_Me_Spring_Lagrangian_Force_Densities(ds,Nb,xLag,yLag,springs):
    ''' Computes the Lagrangian spring force densities
    
    Args:
        ds:
        Nb:
        xLag:
        yLag:
        springs:
        
    Returns:
        fx:
        fy:'''

    #leng = Tx.shape[0]               # # of Lagrangian Pts.

    Nsprings = springs.shape[0]  # # of Springs
    id_Master = springs[:,0].astype('int')   # MASTER NODE Spring Connections
    id_Slave = springs[:,1].astype('int')    # SLAVE NODE Spring Connections
    K_Vec = springs[:,2]  # Stores spring stiffness associated with each spring
    RL_Vec = springs[:,3] # Stores spring resting length associated with each spring
    alpha = springs[:,4]  # Stores deg. of non-linearity of spring

    fx = np.zeros(Nb)                 # Initialize storage for x-forces
    fy = np.zeros(Nb)                 # Initialize storage for y-forces
    
    dx = xLag[id_Slave] - xLag[id_Master] # x-Distance btwn slave and master node
    dy = yLag[id_Slave] - yLag[id_Master] # y-Distance btwn slave and master node
    
    #sF_x = 0.5 * (alpha+1.0) * K_Vec * (np.sqrt(dx**2 + dy**2)-RL_Vec) * (dx/np.sqrt(dx**2+dy**2))
    #sF_y = 0.5 * (alpha+1.0) * K_Vec * (np.sqrt(dx**2 + dy**2)-RL_Vec) * (dy/np.sqrt(dx**2+dy**2))

    sF_x = 0.5 * (alpha+1.0) * K_Vec * np.power( (np.sqrt(dx**2 + dy**2)-RL_Vec), alpha ) * (dx/np.sqrt(dx**2+dy**2))
    sF_y = 0.5 * (alpha+1.0) * K_Vec * np.power( (np.sqrt(dx**2 + dy**2)-RL_Vec), alpha ) * (dy/np.sqrt(dx**2+dy**2))
    #import pdb; pdb.set_trace()
    for ii in range(Nsprings):
        
        fx[id_Master[ii]] += sF_x[ii]  # Sum total forces for node,
                        # i in x-direction (this is MASTER node for this spring)
        fy[id_Master[ii]] += sF_y[ii]  # Sum total forces for node, 
                        # i in y-direction (this is MASTER node for this spring)
        
        fx[id_Slave[ii]] -= sF_x[ii]    # Sum total forces for node, 
                        # i in x-direction (this is SLAVE node for this spring)
        fy[id_Slave[ii]] -= sF_y[ii]    # Sum total forces for node, 
                        # i in y-direction (this is SLAVE node for this spring)

    # MIGHT NOT NEED THESE!
    #fx = fx/ds**2
    #fy = dy/ds**2
    
    return (fx, fy)








###########################################################################
#
# FUNCTION computes the Lagrangian DAMPED SPRING Force Densities .
#
###########################################################################

def give_Me_Damped_Springs_Lagrangian_Force_Densities(ds,Nb,\
            xLag,yLag,d_Springs,xLag_P,yLag_P,dt):
    ''' Computes the Lagrangian spring force densities
    
    Args:
        ds:         Lagrangian spacing
        Nb:         Number of lagrangian points
        xLag:       current xLag positions
        yLag:       current yLag positions
        d_Springs:  holds all information about damped springs
        xLag_P:     previous time-step's xLag positions
        yLag_P:     previous time-step's yLag positions
        dt:         time-step
        
    Returns:
        fx:
        fy:'''

    #leng = Tx.shape[0]               # # of Lagrangian Pts.

    Nsprings = d_Springs.shape[0]              # # of Damped Springs
    id_Master = d_Springs[:,0].astype('int')   # MASTER NODE DAMPED Spring Connections
    id_Slave = d_Springs[:,1].astype('int')    # SLAVE NODE DAMPED Spring Connections
    K_Vec = d_Springs[:,2]  # Stores spring stiffness associated with each Damped spring
    RL_Vec = d_Springs[:,3] # Stores spring resting length associated with each Damped spring
    b_Vec = d_Springs[:,4]  # Store damping coefficient for each Damped spring

    fx = np.zeros(Nb)                 # Initialize storage for x-forces
    fy = np.zeros(Nb)                 # Initialize storage for y-forces
    
    dx = xLag[id_Slave] - xLag[id_Master] # x-Distance btwn slave and master node
    dy = yLag[id_Slave] - yLag[id_Master] # y-Distance btwn slave and master node
    
    dV_x = ( xLag[id_Master] - xLag_P[id_Master] ) / dt  # x-Velocity between current and prev. steps
    dV_y = ( yLag[id_Master] - yLag_P[id_Master] ) / dt  # y-Velocity between current and prev. steps

    sF_x = K_Vec * (np.sqrt(dx**2 + dy**2)-RL_Vec) * (dx/np.sqrt(dx**2+dy**2)) + b_Vec*dV_x 
    sF_y = K_Vec * (np.sqrt(dx**2 + dy**2)-RL_Vec) * (dy/np.sqrt(dx**2+dy**2)) + b_Vec*dV_y
    
    #import pdb; pdb.set_trace()
    for ii in range(Nsprings):
        
        fx[id_Master[ii]] += sF_x[ii]  # Sum total forces for node,
                        # i in x-direction (this is MASTER node for this spring)
        fy[id_Master[ii]] += sF_y[ii]  # Sum total forces for node, 
                        # i in y-direction (this is MASTER node for this spring)
        
        fx[id_Slave[ii]] -= sF_x[ii]    # Sum total forces for node, 
                        # i in x-direction (this is SLAVE node for this spring)
        fy[id_Slave[ii]] -= sF_y[ii]    # Sum total forces for node, 
                        # i in y-direction (this is SLAVE node for this spring)

    # MIGHT NOT NEED THESE!
    #fx = fx/ds**2
    #fy = dy/ds**2
    
    return (fx, fy)










################################################################################
#
# FUNCTION: computes the Mass Lagrangian Force Densities! 
#
################################################################################

def give_Me_Mass_Lagrangian_Force_Densities(ds,xLag,yLag,masses):
    ''' Computes the mass Lagrangian Force Densities
    
    Args:
        ds:
        xLag:
        yLag:
        massess:
        
    Returns:
        fx_mass:
        fy_mass:
        F_mass:'''

    IDs = masses[:,0].astype('int')   # Stores Lag-Pt IDs in col vector
    xPts= masses[:,1]                 # Original x-Values of x-Mass Pts.
    yPts= masses[:,2]                 # Original y-Values of y-Mass Pts.
    kStiffs = masses[:,3]             # Stores "spring" stiffness parameter

    N_masses = masses.shape[0]            # # of target points!

    fx = np.zeros(xLag.size) # Initialize storage for x-force density from TARGET PTS
    fy = np.zeros(xLag.size) # Initialize storage for y-force density from TARGET PTS

    fx[IDs] = fx[IDs] + kStiffs*( xPts - xLag[IDs] )
    fy[IDs] = fy[IDs] + kStiffs*( yPts - yLag[IDs] ) 
       
    fx_mass = fx # alias only
    fy_mass = fy # alias only
    
    F_Mass = np.empty((fx.shape[0],2))
    F_Mass[:,0] = fx  # Store for updating massive boundary pts
    F_Mass[:,1] = fy  # Store for updating massive boundary pts
    
    return (fx_mass, fy_mass, F_Mass)
    
    
    
################################################################################
#
# FUNCTION: computes the Target-Pt Force Densities! 
#
################################################################################

def give_Me_Target_Lagrangian_Force_Densities(ds,xLag,yLag,targets):
    ''' Computes the Target-Pt Densities
    
    Args:
        ds:
        xLag:
        yLag:
        targets:
        
    Returns:
        fx_target:
        fy_target:'''

    IDs = targets[:,0].astype('int')   # Stores Lag-Pt IDs in col vector
    xPts= targets[:,1]                 # Original x-Values of x-Target Pts.
    yPts= targets[:,2]                 # Original y-Values of y-Target Pts.
    kStiffs = targets[:,3]             # Stores Target Stiffnesses 

    N_targets = targets.shape[0]            # # of target points!

    fx = np.zeros(xLag.size)  # Initialize storage for x-force density from TARGET PTS
    fy = np.zeros(xLag.size)  # Initialize storage for y-force density from TARGET PTS

    fx[IDs] = fx[IDs] + kStiffs*( xPts - xLag[IDs] )
    fy[IDs] = fy[IDs] + kStiffs*( yPts - yLag[IDs] )


    fx_target = fx  # alias only
    fy_target = fy  # alias only

    # MIGHT NOT NEED THESE!
    #fx_target = fx/ds**2
    #fy_target = fy/ds**2

    return fx_target, fy_target
    

    
###########################################################################
#
# FUNCTION computes the Lagrangian BEAM Force Densities 
#
###########################################################################

def give_Me_Beam_Lagrangian_Force_Densities(ds,Nb,xLag,yLag,beams):
    '''Computes the Lagrangian BEAM Force Densities
    
    Args:
        ds:
        Nb:
        xLag:
        yLag:
        beams:
        
    Returns:
        fx:
        fy:'''


    Nbeams = beams.shape[0]     # # of Beams
    pts_1 = beams[:,0].astype('int')  # 1ST NODES for BEAM
    pts_2 = beams[:,1].astype('int')  # MIDDLE NODES (2ND Node) for BEAM
    pts_3 = beams[:,2].astype('int')  # 3RD NODES for BEAM
    K_Vec = beams[:,3]  # Stores spring stiffness associated with each spring
    C_Vec = beams[:,4]  # Stores spring resting length associated with each spring

    fx = np.zeros(Nb)                # Initialize storage for x-forces
    fy = np.zeros(Nb)                # Initialize storage for y-forces
    
    Xp = xLag[pts_1]          # xPt of 1ST Node Pt. in beam
    Xq = xLag[pts_2]          # xPt of 2ND (MIDDLE) Node Pt. in beam
    Xr = xLag[pts_3]          # xPt of 3RD Node Pt. in beam
    
    Yp = yLag[pts_1]         # yPt of 1ST Node Pt. in beam
    Yq = yLag[pts_2]         # yPt of 2ND (MIDDLE) Node Pt. in beam
    Yr = yLag[pts_3]         # yPt of 3RD Node Pt. in beam
    
    bF_x =  K_Vec * ( (Xr-Xq)*(Yq-Yp) - (Yr-Yq)*(Xq-Xp) - C_Vec ) * (  (Yq-Yp) + (Yr-Yq) )
    bF_y = -K_Vec * ( (Xr-Xq)*(Yq-Yp) - (Yr-Yq)*(Xq-Xp) - C_Vec ) * (  (Xr-Xq) + (Xq-Xp) )
    
    #bF_x = -K_Vec * ( -(Xr-Xq)*(Yq-Yp) + (Yr-Yq)*(Xq-Xp) - C_Vec ) * (  (Yq-Yp) + (Yr-Yq) )
    #bF_y = -K_Vec * ( -(Xr-Xq)*(Yq-Yp) + (Yr-Yq)*(Xq-Xp) - C_Vec ) * (  -(Xr-Xq) - (Xq-Xp) )
    
    # Want to preserve the ability for the same id_2 to be assigned twice, so
    #   need a loop here.
    for ii in range(Nbeams):
        fx[pts_2[ii]] += bF_x[ii]  # Sum total forces for middle node,
                            # in x-direction (this is MIDDLE node for this beam)
        fy[pts_2[ii]] += bF_y[ii]  # Sum total forces for middle node, 
                            # in y-direction (this is MIDDLE node for this beam)

    # MIGHT NOT NEED THESE!
    #fx = fx/ds**2
    #fy = dy/ds**2

    return (fx, fy)
    
    
    
################################################################################
#
# FUNCTION computes the Delta-Function Approximations 
#
################################################################################

def give_Me_Delta_Function_Approximations_For_Force_Calc(x,y,grid_Info,xLag,yLag):
    '''Computes the Delta-Function Approximations
    
    Args:
        x: np.arange
        y: np.arange
        grid_Info: list
        xLag:
        yLag:
        
    Returns:
        delta_X:
        delta_Y:'''

    # Grid Info
    Nx =   grid_Info['Nx']
    Ny =   grid_Info['Ny']
    Lx =   grid_Info['Lx']
    Ly =   grid_Info['Ly']
    dx =   grid_Info['dx']
    dy =   grid_Info['dy']
    supp = grid_Info['supp']
    Nb =   grid_Info['Nb']

    # Find the indices of the points (xi, yj) where the 1D delta functions are 
    # non-zero in EULERIAN FRAME
    # NOTE: THESE ARE 2D ARRAYS!
    indX = give_1D_NonZero_Delta_Indices(xLag, Nx, dx, supp)
    indY = give_1D_NonZero_Delta_Indices(yLag, Ny, dy, supp) #Python: remove transpose

    # Matrix of possible indices, augmented by "supp"-copies to perform subtractions later in LAGRANGIAN FRAME
    indLagAux = np.arange(Nb)
    # Create copies of indLagAux row vector, stacked vertically.
    # Then transpose to get column vectors
    ind_Lag = np.tile(indLagAux,(supp,1)).T.astype('int')


    # Compute distance between Eulerian Pts. and Lagrangian Pts. by passing correct indices for each
    distX = give_Eulerian_Lagrangian_Distance(x[indX.astype('int')],\
        xLag[ind_Lag],Lx)
    distY = give_Eulerian_Lagrangian_Distance(y[indY.astype('int')],\
        yLag[ind_Lag],Ly) #Python: remove transpose - matches indY
                                        # Result is now tranpose of original code


    # Initialize delta_X and delta_Y matrices for storing delta-function info for each Lag. Pt.
    delta_X = np.zeros((Nb, Nx))
    delta_Y = np.zeros((Ny, Nb))
    

    delta_1D_x = give_Delta_Kernel(distX, dx)
    delta_1D_y = give_Delta_Kernel(distY, dy)

    delta_X[ind_Lag,indX.astype('int')] = delta_1D_x
    delta_Y[indY.astype('int'),ind_Lag] = delta_1D_y
    row,col = ind_Lag.shape

    return (delta_X, delta_Y)
