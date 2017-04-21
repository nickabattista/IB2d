#-------------------------------------------------------------------------------------------------------------------#
#
# IB2d is an Immersed Boundary Code (IB) for solving fully coupled  
#     fluid-structure interaction models. This version of the code is based off of
#    Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
#
# Author: Nicholas A. Battista, Christopher Strickland
# Email:  nick.battista@unc.edu
# Date Created: May 27th, 2015
# Institution: UNC-CH
#
# This code is capable of creating Lagrangian Structures using:
#     1. Springs
#     2. Beams (*torsional springs)
#     3. Target Points
#    4. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-Tension)")
#
# One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting lengths, etc
# 
# There are a number of built in Examples, mostly used for teaching purposes. 
# 
# If you would like us to add a specific muscle model, please let Nick (nick.battista@unc.edu) know.
#
#--------------------------------------------------------------------------------------------------------------------#

#import numpy as np
from read_vtk_data import read_Force_Scalar_Data_From_vtk

#################################################################################
#
# FUNCTION: gives (x,y) positions of the immersed boundary at a single step
#          
################################################################################

def import_Lagrangian_Force_Data(path,numSim):
 
    # read in x directed force #
    strChoice = 'fX_Lag'
    fX_Lag = read_Force_Scalar_Data_From_vtk(path, numSim, strChoice)

    # read in y directed force #
    strChoice = 'fY_Lag'
    fY_Lag = read_Force_Scalar_Data_From_vtk(path, numSim, strChoice)

    # read in y directed force #
    strChoice = 'fMag'
    fLagMag = read_Force_Scalar_Data_From_vtk(path, numSim, strChoice)

    # read in mag. of normal force #
    strChoice = 'fNorm'
    fLagNorm = read_Force_Scalar_Data_From_vtk(path, numSim, strChoice)

    # read in mag. of tangential forces #
    strChoice = 'fTan'
    fLagTan = read_Force_Scalar_Data_From_vtk(path, numSim, strChoice)

    return fX_Lag, fY_Lag, fLagMag, fLagNorm, fLagTan

