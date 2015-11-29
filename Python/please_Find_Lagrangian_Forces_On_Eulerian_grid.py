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
 
################################################################################
#
# FUNCTION: Computes the components of the force term in Navier-Stokes from
#           deformations of the boundary of the immersed boundary
#
################################################################################


def please_Find_Lagrangian_Forces_On_Eulerian_grid(dt, current_time, xLag, yLag,\
    xLag_P,yLag_P, x, y, grid_Info, model_Info, springs, targets, beams,\
    muscles, muscles3, masses):
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
    Nx =    grid_Info[0] # # of Eulerian pts. in x-direction
    Ny =    grid_Info[1] # # of Eulerian pts. in y-direction
    Lx =    grid_Info[2] # Length of Eulerian grid in x-coordinate
    Ly =    grid_Info[3] # Length of Eulerian grid in y-coordinate
    dx =    grid_Info[4] # Spatial-size in x
    dy =    grid_Info[5] # Spatial-size in y
    supp =  grid_Info[6] # Delta-function support
    Nb =    grid_Info[7] # # of Lagrangian pts. 
    ds =    grid_Info[8] # Lagrangian spacing


    # Model Potential Forces #
    springs_Yes = model_Info[0]        # Springs: 0 (for no) or 1 (for yes) 
    target_pts_Yes = model_Info[2]     # Target_Pts: 0 (for no) or 1 (for yes)
    beams_Yes = model_Info[4]          # Beams: 0 (for no) or 1 (for yes)
    muscle_LT_FV_Yes = model_Info[6]   # Length-Tension/Force-Velocity Muscle: 0 (for no) or 1 (for yes) (Length/Tension - Hill Model)
    muscle_3_Hill_Yes = model_Info[7]  # 3-Element Hill Model: 0 (for no) or 1 (for yes) (3 Element Hill + Length-Tension/Force-Velocity)
    mass_Yes = model_Info[10]          # Mass Pts: 0 (for no) or 1 (for yes)


    #
    # Compute MUSCLE LENGTH-TENSION/FORCE-VELOCITY #
    #(if using combined length/tension-Hill model) #
    #
    
    
    if ( muscle_LT_FV_Yes == 1):
        fx_muscles, fy_muscles = give_Muscle_Force_Densities(Nb,xLag,yLag,\
            xLag_P,yLag_P,muscles,current_time,dt)
    else:
        fx_muscles = zeros(xLag.size)
        fy_muscles = zeros(xLag.size)

    

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



    # SUM TOTAL FORCE DENSITY! #
    fx = fx_springs + fx_target + fx_beams + fx_muscles + fx_muscles3 + fx_mass
    fy = fy_springs + fy_target + fy_beams + fy_muscles + fy_muscles3 + fy_mass


    # SAVE LAGRANGIAN FORCES
    F_Lag = np.empty((fx.size,2))
    F_Lag[:,1] = fx
    F_Lag[:,2] = fy
        

    # Give me delta-function approximations!
    delta_X delta_Y = give_Me_Delta_Function_Approximations_For_Force_Calc(x,y,\
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
    m_1 = muscles[:,0]     # Initialize storage for MASTER NODE Muscle Connection
    m_2 = muscles[:,1]     # Initialize storage for SLAVE NODE Muscle Connection
    LFO_Vec = muscles[:,2] # Stores length for max. muscle tension
    SK_Vec = muscles[:,3]  # Stores muscle constant
    a_Vec = muscles[:,4]   # Stores Hill Parameter, a
    b_Vec = muscles[:,5]   # Stores Hill Parameter, b
    FMAX_Vec = muscles[:,6]# Stores Force-Maximum for Muscle

    fx = np.zeros(Nb)                 # Initialize storage for x-forces
    fy = np.zeros(Nb)                 # Initialize storage for y-forces

    for ii in range(Nmuscles):
        
        id_Master = int(m_1[ii])     # Master Node index for i-th muscle
        id_Slave = int(m_2[ii])      # Slave Node index for i-th muscle
        LFO = LFO_Vec[ii]            # Length for max. muscle tension for i-th muscle
        sk = SK_Vec[ii]              # Muscle constant for i-th muscle
        a = a_Vec[ii]                # Hill parameter, a, for i-th muscle
        b = b_Vec[ii]                # Hill parameter, b, for i-th muscle
        Fmax = FMAX_Vec[ii]          # Force-Maximum for i-th muscle
        
        xPt = xLag[id_Master]        # x-Pt of interest at the moment to drive 
                                     #  muscle contraction
        
        dx = xLag[id_Slave] - xLag[id_Master]  # x-Distance btwn slave and master node
        dy = yLag[id_Slave] - yLag[id_Master]  # y-Distance btwn slave and master node
        LF = sqrt( dx**2 + dy**2 )             # Euclidean DISTANCE between master and slave node
        
        
        dx_P = xLag_P[id_Slave] - xLag_P[id_Master]  # x-Distance btwn slave and master node
        dy_P = yLag_P[id_Slave] - yLag_P[id_Master]  # y-Distance btwn slave and master node
        LF_P = sqrt( dx_P**2 + dy_P**2 )             # Euclidean DISTANCE between master and slave node
        
        v =  abs(LF-LF_P)/dt        # How fast the muscle is contracting/expanding
        
        
        # Find actual muscle activation magnitude
        # This function is user defined and included with main2d
        # BE CAREFUL!! xLag is MUTABLE. Pass it by value instead of reference.
        Fm = give_Muscle_Activation(v,LF,LFO,sk,a,b,Fmax,current_time,xPt,\
                                    np.array(xLag))
        
        mF_x = Fm*(dx/LF)           # cos(theta) = dx / LF;
        mF_y = Fm*(dy/LF)           # sin(theta) = dy / LF;
        
        fx[id_Master] = fx[id_Master] + mF_x  # Sum total forces for node,
                        # i in x-direction (this is MASTER node for this spring)
        fy[id_Master] = fy[id_Master] + mF_y  # Sum total forces for node,
                        # i in y-direction (this is MASTER node for this spring)
        
        fx[id_Slave] = fx[id_Slave] - mF_x    # Sum total forces for node,
                        # i in x-direction (this is SLAVE node for this spring)
        fy[id_Slave] = fy[id_Slave] - mF_y    # Sum total forces for node,
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
    m_1 = muscles3[:,0]      # Initialize storage for MASTER NODE Muscle Connection
    m_2 = muscles3[:,1]      # Initialize storage for SLAVE NODE Muscle Connection
    LFO_Vec = muscles3[:,2]  # Stores length for max. muscle tension
    SK_Vec = muscles3[:,3]   # Stores muscle constant
    a_Vec = muscles3[:,4]    # Stores Hill Parameter, a
    b_Vec = muscles3[:,5]    # Stores Hill Parameter, b
    FMAX_Vec = muscles3[:,6] # Stores Force-Maximum for Muscle

    fx = np.zeros(Nb)                 # Initialize storage for x-forces
    fy = np.zeros(Nb)                 # Initialize storage for y-forces

    for ii in range(Nmuscles):
        
        id_Master = m_1[ii]          # Master Node index for i-th muscle
        id_Slave = m_2[ii]           # Slave Node index for i-th muscle
        LFO = LFO_Vec[ii]            # Length for max. muscle tension for i-th muscle
        sk = SK_Vec[ii]              # Muscle constant for i-th muscle
        a = a_Vec[ii]                # Hill parameter, a, for i-th muscle
        b = b_Vec[ii]                # Hill parameter, b, for i-th muscle
        Fmax = FMAX_Vec[ii]          # Force-Maximum for i-th muscle
        
        xPt = xLag[ id_Master ]     # x-Pt of interest at the moment to drive muscle contraction
        
        dx = xLag[id_Slave] - xLag[id_Master] # x-Distance btwn slave and master node
        dy = yLag[id_Slave] - yLag[id_Master] # y-Distance btwn slave and master node
        LF = sqrt( dx**2 + dy**2 )            # Euclidean DISTANCE between 
                                              #   master and slave node
        
        
        dx_P = xLag_P[id_Slave] - xLag_P[id_Master] # x-Distance btwn slave and master node
        dy_P = yLag_P[id_Slave] - yLag_P[id_Master] # y-Distance btwn slave and master node
        LF_P = sqrt( dx_P**2 + dy_P**2 )            # Euclidean DISTANCE between 
                                                    #   master and slave node
        
        v =  abs(LF-LF_P)/dt         # How fast the muscle is contracting/expanding 

        # Find actual muscle activation magnitude
        # This function is user defined and included with main2d
        # BE CAREFUL!! xLag is MUTABLE. Pass it by value instead of reference.
        Fm = give_3_Element_Muscle_Activation(v,LF,LFO,sk,a,b,Fmax,current_time,\
            xPt,np.array(xLag))
        
        mF_x = Fm*(dx/LF)           # cos(theta) = dx / LF;
        mF_y = Fm*(dy/LF)           # sin(theta) = dy / LF;
        
        fx[id_Master] = fx[id_Master] + mF_x  # Sum total forces for node,
                        # i in x-direction (this is MASTER node for this spring)
        fy[id_Master] = fy[id_Master] + mF_y  # Sum total forces for node, 
                        # i in y-direction (this is MASTER node for this spring)
        
        fx[id_Slave] = fx[id_Slave] - mF_x    # Sum total forces for node, 
                        # i in x-direction (this is SLAVE node for this spring)
        fy[id_Slave] = fy[id_Slave] - mF_y    # Sum total forces for node, 
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
    sp_1 = springs[:,0]   # Initialize storage for MASTER NODE Spring Connection
    sp_2 = springs[:,1]   # Initialize storage for SLAVE NODE Spring Connection
    K_Vec = springs[:,2]  # Stores spring stiffness associated with each spring
    RL_Vec = springs[:,3] # Stores spring resting length associated with each spring

    fx = np.zeros(Nb)                 # Initialize storage for x-forces
    fy = np.zeros(Nb)                 # Initialize storage for y-forces

    for ii in range(Nsprings):
        
        id_Master = sp_1[ii]          # Master Node index
        id_Slave = sp_2[ii]           # Slave Node index
        k_Spring = K_Vec[ii]          # Spring stiffness of i-th spring
        L_r = RL_Vec[ii]              # Resting length of i-th spring
        
        dx = xLag[id_Slave] - xLag[id_Master] # x-Distance btwn slave and master node
        dy = yLag[id_Slave] - yLag[id_Master] # y-Distance btwn slave and master node
        
        sF_x = k_Spring * ( sqrt(dx**2 + dy**2) - L_r ) * ( dx / sqrt(dx**2+dy**2) )
        sF_y = k_Spring * ( sqrt(dx**2 + dy**2) - L_r ) * ( dy / sqrt(dx**2+dy**2) )
        
        fx[id_Master] = fx[id_Master] + sF_x  # Sum total forces for node,
                        # i in x-direction (this is MASTER node for this spring)
        fy[id_Master] = fy[id_Master] + sF_y  # Sum total forces for node, 
                        # i in y-direction (this is MASTER node for this spring)
        
        fx[id_Slave] = fx[id_Slave] - sF_x    # Sum total forces for node, 
                        # i in x-direction (this is SLAVE node for this spring)
        fy[id_Slave] = fy[id_Slave] - sF_y    # Sum total forces for node, 
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

    IDs = masses[:,0]                 # Stores Lag-Pt IDs in col vector
    xPts= masses[:,1]                 # Original x-Values of x-Mass Pts.
    yPts= masses[:,2]                 # Original y-Values of y-Mass Pts.
    kStiffs = masses[:,3]             # Stores "spring" stiffness parameter

    N_masses = masses.shape[0]            # # of target points!

    fx = np.zeros(xLag.size) # Initialize storage for x-force density from TARGET PTS
    fy = np.zeros(xLag.size) # Initialize storage for y-force density from TARGET PTS

    for ii in range(N_masses):
        fx[int(IDs[ii])] = fx[int(IDs[ii])] + kStiffs[ii]*( xPts[ii] - xLag[int(IDs[ii])] )
        fy[int(IDs[ii])] = fy[int(IDs[ii])] + kStiffs[ii]*( yPts[ii] - yLag[int(IDs[ii])] ) 
       
    fx_mass = fx # alias only
    fy_mass = fy # alias only
    
    F_Mass = np.empty((fx.shape[0],2))
    F_Mass[:,1] = fx  # Store for updating massive boundary pts
    F_Mass[:,2] = fy  # Store for updating massive boundary pts
    
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

    IDs = targets[:,0]                 # Stores Lag-Pt IDs in col vector
    xPts= targets[:,1]                 # Original x-Values of x-Target Pts.
    yPts= targets[:,2]                 # Original y-Values of y-Target Pts.
    kStiffs = targets[:,3]             # Stores Target Stiffnesses 

    N_targets = targets.shape[0]            # # of target points!

    fx = np.zeros(xLag.size)  # Initialize storage for x-force density from TARGET PTS
    fy = np.zeros(xLag.size)  # Initialize storage for y-force density from TARGET PTS

    for ii in range(N_targets):
        fx[int(IDs[ii])] = fx[int(IDs[ii])] + kStiffs[ii]*( xPts[ii] - xLag[int(IDs[ii])] )
        fy[int(IDs[ii])] = fy[int(IDs[ii])] + kStiffs[ii]*( yPts[ii] - yLag[int(IDs[ii])] ) 


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
    pts_1 = beams[:,0]  # Initialize storage for 1ST NODE for BEAM
    pts_2 = beams[:,1]  # Initialize storage for MIDDLE NODE (2ND Node) for BEAM
    pts_3 = beams[:,2]  # Initialize storage for 3RD NODE for BEAM
    K_Vec = beams[:,3]  # Stores spring stiffness associated with each spring
    C_Vec = beams[:,4]  # Stores spring resting length associated with each spring

    fx = np.zeros(Nb)                # Initialize storage for x-forces
    fy = np.zeros(Nb)                # Initialize storage for y-forces

    for ii in range(Nbeams):
        
        id_1 = int(pts_1[ii])  # 1ST Node index
        id_2 = int(pts_2[ii])  # (MIDDLE) 2nd Node index -> 
                               #           index that gets force applied to it!
        id_3 = int(pts_3[ii])  # 3RD Node index
        k_Beam = K_Vec[ii]     # Beam stiffness of i-th spring
        C = C_Vec[ii]          # Curvature of the beam between these three nodes 
        
        Xp = xLag[id_1]          # xPt of 1ST Node Pt. in beam
        Xq = xLag[id_2]          # xPt of 2ND (MIDDLE) Node Pt. in beam
        Xr = xLag[id_3]          # xPt of 3RD Node Pt. in beam
        
        Yp = yLag[id_1]         # yPt of 1ST Node Pt. in beam
        Yq = yLag[id_2]         # yPt of 2ND (MIDDLE) Node Pt. in beam
        Yr = yLag[id_3]         # yPt of 3RD Node Pt. in beam
        
        bF_x =  k_Beam * ( (Xr-Xq)*(Yq-Yp) - (Yr-Yq)*(Xq-Xp) - C ) * (  (Yq-Yp) + (Yr-Yq) )
        bF_y = -k_Beam * ( (Xr-Xq)*(Yq-Yp) - (Yr-Yq)*(Xq-Xp) - C ) * (  (Xr-Xq) + (Xq-Xp) )
        
        #bF_x = -k_Beam * ( -(Xr-Xq)*(Yq-Yp) + (Yr-Yq)*(Xq-Xp) - C ) * (  (Yq-Yp) + (Yr-Yq) )
        #bF_y = -k_Beam * ( -(Xr-Xq)*(Yq-Yp) + (Yr-Yq)*(Xq-Xp) - C ) * (  -(Xr-Xq) - (Xq-Xp) )
        
        fx[id_2] = fx[id_2] + bF_x  # Sum total forces for middle node,
                            # in x-direction (this is MIDDLE node for this beam)
        fy[id_2] = fy[id_2] + bF_y  # Sum total forces for middle node, 
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
    Nx =   grid_Info[0]
    Ny =   grid_Info[1]
    Lx =   grid_Info[2]
    Ly =   grid_Info[3]
    dx =   grid_Info[4]
    dy =   grid_Info[5]
    supp = grid_Info[6]
    Nb =   grid_Info[7]

    # Find the indices of the points (xi, yj) where the 1D delta functions are 
    # non-zero in EULERIAN FRAME
    # NOTE: THESE ARE 2D ARRAYS!
    indX = give_1D_NonZero_Delta_Indices(xLag, Nx, dx, supp)
    indY = give_1D_NonZero_Delta_Indices(yLag, Ny, dy, supp) #Python: remove transpose

    # Matrix of possible indices, augmented by "supp"-copies to perform subtractions later in LAGRANGIAN FRAME
    indLagAux = np.arange(Nb)
    # Create copies of indLagAux row vector, stacked vertically.
    # Then transpose to agree with current column major code
    ind_Lag = np.vstack(indLagAux for ii in range(supp)).T


    # Compute distance between Eulerian Pts. and Lagrangian Pts. by passing correct indices for each
    distX = give_Eulerian_Lagrangian_Distance(x[indX.astype('int')],\
        xLag[ind_Lag.astype('int')],Lx)
    distY = give_Eulerian_Lagrangian_Distance(y[indY.astype('int')],\
        yLag[ind_Lag.astype('int')],Ly) #Python: remove transpose - matches indY
                                        # Result is now tranpose of original code


    # Initialize delta_X and delta_Y matrices for storing delta-function info for each Lag. Pt.
    delta_X = np.zeros((Nb, Nx))
    delta_Y = np.zeros((Ny, Nb))
    

    delta_1D_x = give_Delta_Kernel(distX, dx)
    delta_1D_y = give_Delta_Kernel(distY, dy)


    row,col = ind_Lag.shape
    for ii in range(row):
        for jj in range(col):
            
            # Get Eulerian/Lagrangian indices to use for saving non-zero delta-function values
            xID = indX[ii,jj]
            indy= ind_Lag[ii,jj]
            yID = indY[ii,jj] #Python: now indicies don't need to be switched
            
            # Store non-zero delta-function values into delta_X / delta_Y matrices at correct indices!
            delta_X[indy,xID] = delta_1D_x[ii,jj]
            delta_Y[yID,indy] = delta_1D_y[ii,jj] #now return transpose

    return (delta_X delta_Y)
    
    
    
###########################################################################
#
# FUNCTION finds the indices on the Eulerian grid where the 1D Dirac-delta
# kernel is possibly non-zero is x-dimension.
#
###########################################################################

def give_1D_NonZero_Delta_Indices(lagPts_j, N, dx, supp):
    '''Find the indices on Eulerian grid where 1D delta kernel is poss. non-zero
    
    Args:
        lagPts_j: 2D array of lagrangian pts for specific coordinate, j= x or y.
        N: # spatial resolution of Eulerian grid in each dimension
        dx: Spatial step-size on Eulerian (fluid) grid
        supp: Size of support of the Dirac-delta kernel (should be even)
        
    Returns:
        indices:''' 

    # Finds the index of the lower left Eulerian pt. to Lagrangian pt..
    ind_Aux = np.floor(lagPts_j/dx + 1)

    # Get all the different x indices that must be considered.
    indices = np.vstack(ind_Aux for ii in range(supp)) #ind_Aux is each row
    #
    for ii in range(supp):
        indices[ii,:] = indices[ii,:] + -supp/2+1+ii

    # Translate indices between {0,2,..,N-1}
    indices = indices-1 % N
    
    # Transpose to work with current column-major code
    return indices.T
    

    
################################################################################
#
# FUNCTION distance between Eulerian grid data, x, and Lagrangian grid data, y, 
#          at specifed pts typically and makes sure the distance are [0,L] accordingly.
#
################################################################################

def give_Eulerian_Lagrangian_Distance(x, y, L):
    ''' Distance between Eulerian grid data x and Lagrangian grid data y within [0,L]
    
    Args:
        x: 2D ndarray (Eulerian data)
        y: 2D ndarray (Lagrangian data)
        L: Length of domain, [0,L]
        
    Returns:
        distance:'''

    row,col = x.shape
    distance = np.abs( x - y )
    for ii in range(row):
        for jj in range(col):
            distance[ii,jj] = min( distance(ii,jj), L-distance(ii,jj) ) 
            #Note: need to make sure that taking correct value
            #CHRISTOPHER'S NOTE: This looks like an error to me!

    return distance

    

    
###########################################################################
#
# FUNCTION: computes a discrete approx. to a 1D Dirac-delta function over a
# specified matrix, x, and spatial step-size, dx. It will have support in
# [x-2dx, x+2dx]
#
###########################################################################

def give_Delta_Kernel(x,dx):
    ''' Compute discrete approx. to 1D delta function over x, dx.
    Support is in [x-2dx, x+2dx].
    
    Args:
        x:  Values in which the delta function will be evaulated
        dx: Spatial step-size of grid
        
    Returns:
        delta:'''

    # Computes Dirac-delta Approximation.
    RMAT = np.abs(x)/dx

    #Initialize delta
    delta = np.array(RMAT)

    #Loops over to find delta approximation
    row,col = x.shape
    for ii in range(row):
        for jj in range(col):
            
            r = RMAT[ii,jj]
            
            if r<1:
                delta[ii,jj] = ( (3 - 2*r + sqrt(1 + 4*r - 4*r*r) ) / (8*dx) )
            elif ( (r<2) and (r>=1) ):
                delta[ii,jj] = ( (5 - 2*r - sqrt(-7 + 12*r - 4*r*r) ) / (8*dx) )

    return delta


