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

import numpy as np
import os

#################################################################################
#
# FUNCTION: Reads in (x,y) positions of the immersed boundary from .vtk format
#          
################################################################################

def read_Eulerian_Velocity_Field_vtk(path,simNums):
    
    analysis_path = os.getcwd();  # Stores current directory path        
           
    os.chdir(path);               # cd's into viz_IB2d folder
    
    filename = 'u.' + str(simNums) + '.vtk';

    # Stores grid resolution from .vtk file
    Nx = np.genfromtxt(filename, skip_header=5, usecols=(1),max_rows=1);    
    
    # Stores desired Eulerian data
    e_data_X = np.genfromtxt(filename, skip_header=13, usecols=range(0,3*int(Nx),3),max_rows=int(Nx))
    e_data_Y = np.genfromtxt(filename, skip_header=13, usecols=range(1,3*int(Nx),3),max_rows=int(Nx))

    os.chdir(analysis_path);     # Path to working directory

    return e_data_X, e_data_Y