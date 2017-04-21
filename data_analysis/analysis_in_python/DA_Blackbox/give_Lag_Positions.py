#-------------------------------------------------------------------------------------------------------------------#
#
# IB2d is an Immersed Boundary Code (IB) for solving fully coupled  
# 	fluid-structure interaction models. This version of the code is based off of
#	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
#
# Author: Nicholas A. Battista, Christopher Strickland
# Email:  nick.battista@unc.edu
# Date Created: May 27th, 2015
# Institution: UNC-CH
#
# This code is capable of creating Lagrangian Structures using:
# 	1. Springs
# 	2. Beams (*torsional springs)
# 	3. Target Points
#	4. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-Tension)")
#
# One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting lengths, etc
# 
# There are a number of built in Examples, mostly used for teaching purposes. 
# 
# If you would like us to add a specific muscle model, please let Nick (nick.battista@unc.edu) know.
#
#--------------------------------------------------------------------------------------------------------------------#

import numpy as np
from pathlib import Path
from read_vtk_data import read_vtk_Unstructured_Grid_Points

#################################################################################
#
# FUNCTION: gives (x,y) positions of the immersed boundary at a single step
#          
################################################################################

def give_Lag_Positions(path,numSim):
    
    ''' path: (string) gives path to viz_IB2d directory
        numSim: number of simulation output to tac onto path
        
        RETURNS:
            xLag: x-Lagrangian point positions
            yLag: y-Lagrangian point positions'''
          
    filename = Path(path) / ('lagsPts.' + str(numSim) + '.vtk')
                 
    lag_data = read_vtk_Unstructured_Grid_Points(str(filename))

    xLag = lag_data[:,0]
    yLag = lag_data[:,1]

    return xLag, yLag