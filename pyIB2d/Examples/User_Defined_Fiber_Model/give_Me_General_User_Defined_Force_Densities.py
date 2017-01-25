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
from math import cos, sin, pi, sqrt

######################################################################################
#
# FUNCTION: Models the force interaction that is user-defined.
#
#           Note:
#                   1. User has complete control for artibtrary input file
#                      format and number of parameters
#                   2. Those parameters get passed to the matrix called
#                      "general_force"
#                   3. User gets options of including other things into
#                      force function outside of inputted parameters, e.g.,
#                      current_time, dt, previous location of Lag. Pts, etc.
#
######################################################################################

def give_Me_General_User_Defined_Force_Densities(ds,Nb,xLag,yLag,\
            xLag_P,yLag_P,dt,current_time,general_force):


    #
    # INPUTS:
    #         ds: Lagragian spacing (defined by ds = 0.5*dx)
    #         Nb: # of Lagrangian Pts.
    #         xLag: current x-Lagrangian coordinate positions
    #         yLag: current y-Lagrangian coordinate positions
    #         xLag_P: previous x-Lagrangian coordinate positions
    #         yLag_P: previous y-Lagrangian coordinate positions
    #         dt: time-step value
    #         current_time: current time in simulation
    #         general_force: matrix containing all data from input file

    
    #
    # NOTE: THIS EXAMPLE CREATES A USER-DEFINED FORCE THAT IS JUST AN
    # ORDINARY LINEAR SPRING
    #

    Nsprings = general_force.shape[1]  # # of Force Connections (e.g., springs here)

    id_Master = general_force[0,:].astype('int')   # MASTER NODE Spring Connections
    id_Slave = general_force[1,:].astype('int')    # SLAVE NODE Spring Connections
    K_Vec = general_force[2,:]  # Stores spring stiffness associated with each spring
    RL_Vec = general_force[3,:] # Stores spring resting length associated with each spring

    fx = np.zeros(Nb)                 # Initialize storage for x-forces
    fy = np.zeros(Nb)                 # Initialize storage for y-forces
    
    dx = xLag[id_Slave] - xLag[id_Master] # x-Distance btwn slave and master node
    dy = yLag[id_Slave] - yLag[id_Master] # y-Distance btwn slave and master node

    sF_x = K_Vec * (np.sqrt(dx**2 + dy**2)-RL_Vec)  * (dx/np.sqrt(dx**2+dy**2))
    sF_y = K_Vec * (np.sqrt(dx**2 + dy**2)-RL_Vec)  * (dy/np.sqrt(dx**2+dy**2))

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
    
    return (fx, fy)
