#-------------------------------------------------------------------------------------------------------------------#
#
# IB2d is an Immersed Boundary Code (IB) for solving fully coupled  
# 	fluid-structure interaction models. This version of the code is based off of
#	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
#
# Author: Nicholas A. Battista
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

#import numpy as np
#import sys
import os
from read_Lagrangian_Data_From_vtk import read_Lagrangian_Data_From_vtk
#from math import cos, sin, pi, sqrt
#import matplotlib.pyplot as plt

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
          
    analysis_path = os.getcwd()
                 
    xLag,yLag = read_Lagrangian_Data_From_vtk(path,numSim)

    os.chdir(analysis_path)

    return xLag,yLag