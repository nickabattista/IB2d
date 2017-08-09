#*****************************************************************************#
#***********************************# IB2d #**********************************#
#*****************************************************************************#

# IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
#  	fluid-structure interaction models. This version of the code is based 
# 	off of Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
#
# Author: Nicholas A. Battista
# Email:  nickabattista@gmail.com
# Date Created: May 27th, 2015
# Institution: University of North Carolina at Chapel Hill
# Website: http://battista.web.unc.edu
# GitHub: http://www.github.com/nickabattista
# 
# This code is capable of creating Lagrangian Structures using:
#  	1. Springs
#  	2. Beams (*torsional springs)
#  	3. Target Points
# 	4. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-
#                         Tension)")
#   . Mass Points
# 
# One is able to update those Lagrangian Structure Parameters, e.g., spring constants, resting lengths, etc
#  
# There are a number of built in Examples, mostly used for teaching purposes. 
# 
# If you would like us to add a specific muscle model, please contact Nick (nickabattista@gmail.com) 
# 
# If you use this code for the purposes of teaching, research, or recreation please let Nick know as well :)
    
import numpy as np

###########################################################################
#
# FUNCTION: Reads in Eulerian nodes for applying Boussinesq Approximation
#
###########################################################################    

def get_Inds_Please():
       
    #Reads in the number Boussinesq indices of Eulerian mesh
    
    #Args: 
        
    #Returns:
    #    N: number of Lagrangian points
    #    xLag: x-values of Lagrangian mesh
    #    yLab: y-values of Lagrangian mesh

    filename = 'boussinesq.mesh'  # Name of file to read in
    with open(filename) as f:
        # First line in the file contains the number of Lagrangian points
        # N = int(f.readline().strip())
        # Read in the Lagrangian mesh points
        xInds,yInds = np.loadtxt(f,unpack=True)
    
    xInds = xInds.astype(int)    
    yInds = yInds.astype(int)
    return (xInds,yInds)


############################################################################################
#
# FUNCTION: -computes Boussinesq Forcing Terms to pass into Fluid Solve
#           -user has option to customize
#
############################################################################################

def please_Form_Boussinesq_Forcing_Terms(exp_Coeff,Nx,Ny,gravity_Info):

    # INPUTS:
    # exp_Coeff:    coefficient of expansion, e.g., thermal expansion
    # Nx,Ny:        grid resolution in x,y directions
    # gravity_Info: 
    #               col 0: flag if considering gravity
    #               col 1: x-component of gravity vector (normalized)
    #               col 2: y-component of gravity vector (normalized)

    g = 0.981       # gravitational constant
    mat = np.zeros((Ny,Nx))
    
    xInds,yInds = get_Inds_Please()
    for i in range (0,xInds.shape[0]):
        xInd = xInds[i]
        yInd = yInds[i]
        mat[yInd,xInd] = 1

    fBouss_X = exp_Coeff*g*gravity_Info[1]*mat
    fBouss_Y = exp_Coeff*g*gravity_Info[2]*mat

    return (fBouss_X, fBouss_Y)
    


