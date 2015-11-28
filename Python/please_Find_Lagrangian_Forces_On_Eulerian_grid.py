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
        fx_muscles = zeros(len(xLag))
        fy_muscles = zeros(len(xLag))

    pass

    #
    # Compute 3-ELEMENT HILL MUSCLE MODEL FORCE DENSITIES (if using combined 3-HILL + LT/FV) #
    #
    if ( muscle_3_Hill_Yes == 1)
        [fx_muscles3, fy_muscles3] = give_3_Element_Muscle_Force_Densities(Nb,xLag,yLag,xLag_P,yLag_P,muscles3,current_time,dt);
    else
        fx_muscles3 = zeros(length(xLag),1);
        fy_muscles3 = fx_muscles3;
    end




    # Compute SPRING FORCE DENSITIES (if there are springs!)
    if ( springs_Yes == 1 )
        
        # Compute "Connections Matrix" for what springs are attached to whom     #
        # (and "Connections Stiffness Matrix" and "Connections Restling Lengths" #
        #connects = give_Me_Spring_Connections_Matrix(Nb,Nsprings,sp_1,sp_2,K_Vec,L_Vec); 
        
        # Compute distances between Lag-Pts w/ Springs for Spring-Tension Calc.
        #dLag_x = give_Spring_Lagrangian_Distance(xLag, Lx, springs);
        #dLag_y = give_Spring_Lagrangian_Distance(yLag, Ly, springs);

        # Compute the Lagrangian SPRING tensions!
        #[Tx Ty] = give_Me_Spring_Lagrangian_Tension(Nb,dLag_x,dLag_y,springs);

        # Compute the Lagrangian SPRING force densities!
        [fx_springs, fy_springs] = give_Me_Spring_Lagrangian_Force_Densities(ds,Nb,xLag,yLag,springs);
        
    else
        fx_springs = zeros(Nb,1); #No x-forces coming from springs
        fy_springs = fx_springs;  #No y-forces coming from springs
    end



    # Compute MASS PT FORCE DENSITIES (if there are mass points!)
    if ( mass_Yes == 1)
        # Compute the Lagrangian MASS PT force densities!
        [fx_mass, fy_mass, F_Mass] = give_Me_Mass_Lagrangian_Force_Densities(ds,xLag,yLag,masses); 
    else
        fx_mass = zeros(Nb,1); #No x-forces coming from mass points
        fy_mass = fx_mass;     #No y-forces coming from mass points
        F_Mass = 0;            #Dummy to pass along  
    end



    # Compute TARGET FORCE DENSITIES (if there are target points!)
    if ( target_pts_Yes == 1)
        # Compute the Lagrangian TARGET force densities!
        [fx_target, fy_target] = give_Me_Target_Lagrangian_Force_Densities(ds,xLag,yLag,targets); 
        
    else
        fx_target = zeros(Nb,1); #No x-forces coming from target points
        fy_target = fx_target;   #No y-forces coming from target points
    end




    # Compute BEAM FORCE DENSITIES (if there are beams!)
    if ( beams_Yes == 1 )

        # Compute the Lagrangian SPRING force densities!
        [fx_beams, fy_beams] = give_Me_Beam_Lagrangian_Force_Densities(ds,Nb,xLag,yLag,beams);
        
    else
        fx_beams = zeros(Nb,1); #No x-forces coming from beams
        fy_beams = fx_beams;    #No y-forces coming from beams
    end



    # SUM TOTAL FORCE DENSITY! #
    fx = fx_springs + fx_target + fx_beams + fx_muscles + fx_muscles3 + fx_mass;
    fy = fy_springs + fy_target + fy_beams + fy_muscles + fy_muscles3 + fy_mass;


    # SAVE LAGRANGIAN FORCES
    F_Lag(:,1) = fx;
    F_Lag(:,2) = fy;
        

    # Give me delta-function approximations!
    [delta_X delta_Y] = give_Me_Delta_Function_Approximations_For_Force_Calc(x,y,grid_Info,xLag,yLag);


    # Transform the force density vectors into diagonal matrices
    fxds = zeros(Nb,Nb); fyds = zeros(Nb,Nb);
    for i=1:Nb
       fxds(i,i) = fx(i,1)*ds; 
       fyds(i,i) = fy(i,1)*ds;
    end


    # Find Eulerian forces on grids by approximating the line integral, 
    #       F(x,y) = int{ f(s) delta(x - xLag(s)) delta(y - yLag(s)) ds }
    Fx = delta_Y * fxds * delta_X;
    Fy = delta_Y * fyds * delta_X;

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
    # Need np.array() to get assignment by value
    Nmuscles = muscles.shape[0]      # # of Muscles
    m_1 = np.array(muscles[:,0])     # Initialize storage for MASTER NODE Muscle Connection
    m_2 = np.array(muscles[:,1])     # Initialize storage for SLAVE NODE Muscle Connection
    LFO_Vec = np.array(muscles[:,2]) # Stores length for max. muscle tension
    SK_Vec = np.array(muscles[:,3])  # Stores muscle constant
    a_Vec = np.array(muscles[:,4])   # Stores Hill Parameter, a
    b_Vec = np.array(muscles[:,5])   # Stores Hill Parameter, b
    FMAX_Vec = np.array(muscles[:,6])# Stores Force-Maximum for Muscle

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