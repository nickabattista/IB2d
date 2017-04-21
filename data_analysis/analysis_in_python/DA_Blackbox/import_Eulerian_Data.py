#-------------------------------------------------------------------------------------------------------------------#
#
# IB2d is an Immersed Boundary Code (IB) for solving fully coupled  
#     fluid-structure interaction models. This version of the code is based off of
#    Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
#
# Author: Nicholas A. Battista
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
#import sys
#import os
from read_vtk_data import read_Eulerian_Data_From_vtk

#################################################################################
#
# FUNCTION: gives (x,y) positions of the immersed boundary at a single step
#          
################################################################################

def import_Eulerian_Data(path,numSim):
 
    # read in Vorticity #
    strChoice = 'Omega'; getxy = True
    Omega,x,y = read_Eulerian_Data_From_vtk(path,numSim,strChoice,getxy)

    # read in Pressure #
    strChoice = 'P'; getxy = False
    P = read_Eulerian_Data_From_vtk(path,numSim,strChoice,getxy)

    # read in Velocity Magnitude #
    strChoice = 'uMag'; getxy = False
    uMag = read_Eulerian_Data_From_vtk(path,numSim,strChoice,getxy)
 
    # read in x-directed Velocity Magnitude #
    strChoice = 'uX'; getxy = False
    uX = read_Eulerian_Data_From_vtk(path,numSim,strChoice,getxy)
 
    # read in y-directed Velocity Magnitude #
    strChoice = 'uY'; getxy = False
    uY = read_Eulerian_Data_From_vtk(path,numSim,strChoice,getxy)

    # read in x-directed Forces #
    strChoice = 'Fx'; getxy = False
    Fx = read_Eulerian_Data_From_vtk(path,numSim,strChoice,getxy)

    # read in y-directed Forces #
    strChoice = 'Fy'; getxy = False
    Fy = read_Eulerian_Data_From_vtk(path,numSim,strChoice,getxy)
 
    # read in Velocity Field #
    strChoice = 'u'; getxy = False
    U,V = read_Eulerian_Data_From_vtk(path,numSim,strChoice,getxy)

    return x,y,Omega,P,uMag,uX,uY,U,V,Fx,Fy

