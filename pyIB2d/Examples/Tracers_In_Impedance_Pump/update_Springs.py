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
from math import pi
 
#####################################################################
#
# FUNCTION: updates the resting lengths and/or spring stiffness
#
#####################################################################

def update_Springs(dt,current_time,xLag,yLag,springs):
    ''' Updates the resting lengths and/or spring stiffness
    
    Args:
        dt: time-step
        current_time: current time of simulation
        xLag: Vector of all xLagrangian Pts        
        springs: col 1: starting spring pt (by lag. discretization)
                 col 2: ending spring pt. (by lag. discretization)
                 col 3: spring stiffness
                 col 4: spring resting lengths
                 
    Returns:
        springs:'''

    s_1 = springs[:,0]    # MASTER Node Indices
    s_2 = springs[:,1]    # SLAVE Node Indices
    k_Vec = springs[:,2]  # Springs Stiffnesses Vector
    RL_Vec= springs[:,3]  # Resting Lengths Vectors

    N = xLag.size                   # Gives total number of Lagrangian pts!
    N_springs = springs.shape[0]    # Gives total number of springs!

    freq = 10  # Frequency for pumping
    d = 1.0    # Diameter of Heart Tube
    
    RL_Vec[N-3+10:N-2+20] = d - np.abs( 0.9*d*np.sin( 2*pi*freq*current_time ) )

    springs[:,3] = RL_Vec # Update the springs_info

    return springs

