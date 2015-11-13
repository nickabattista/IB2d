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
import sys
# Path Reference to where Driving code is found #
sys.path.append('../')
import IBM_Driver as Driver

###############################################################################
#
# Function to read in input files from "input2d" -> renames all quantities
# appropriately, just as they are in the input file
#
###############################################################################

def give_Me_input2d_Parameters():
    ''' Function to read in input files from "input2d"
    
    Returns:
        params: ndarray of parameter values
        struct_name: name of structure'''
        
    filename = 'input2d' #Name of file to read in
    
    #This is more sophisticated than what was here before the port.
    #If desired, names could be double checked for consistancy.
    names,params = np.loadtxt(filename,dtype={'names': ('param','value'),\
            'formats': ('S25','f8')},comments=('%','string_name'),\
            delimiter='=',unpack=True)
            
    #now get string_name
    with open('input2d','r') as f:
        for line in f:
            #look for 'string_name'. make sure it comes before a comment char
            if ('string_name' in line) and \
                line.find('string_name') < line.find('%'):
                #split the line by whitespace
                words = line.split()
                #look for the equals sign
                for n,word in enumerate(words):
                    if word == '=':
                        struct_name = words[n+1]
                        break
                break
    
    return (params,struct_name)

###############################################################################
#
# FUNCTION: main2d is the function that gets called to run the code. It
#           itself reads in paramters from the input2d file, and passes
#           them to the IBM_Driver function to run the simulation
#
###############################################################################

def main2d():
    
    '''This is the "main" function, which ets called to run the 
    Immersed Boundary Simulation. It reads in all the parameters from 
    "input2d", and sends them off to the "IBM_Driver" function to actually 
    perform the simulation.'''
    
    # READ-IN INPUT PARAMTERS #
    params,struct_name = give_Me_input2d_Parameters()
    
    # FLUID PARAMETER VALUES STORED #
    mu = params[0]      # Dynamic Viscosity
    rho = params[1]     # Density

    # TEMPORAL INFORMATION VALUES STORED #
    T_final = params[2] # Final simulation time
    dt = params[3]      # Time-step

    # GRID INFO STORED #
    grid_Info = list(range(14))
    grid_Info[0] = params[4]           # num of Eulerian Pts in x-Direction
    grid_Info[1] = params[5]           # num of Eulerian Pts in y-Direction 
    grid_Info[2] = params[6]           # Length of Eulerian domain in x-Direction
    grid_Info[3] = params[7]           # Length of Eulerian domain in y-Direction
    grid_Info[4] = params[6]/params[4] # Spatial step-size in x
    grid_Info[5] = params[7]/params[5] # Spatial step-size in y
    grid_Info[6] = params[8]           # num of pts used in delta-function support 
                                       #    (supp/2 in each direction)
    grid_Info[7] = params[25]          # Print Dump (How often to plot)
    grid_Info[8] = params[26]          # Plot in Matlab? (1=YES,0=NO) 
    grid_Info[9] = params[27]          # Plot LAGRANGIAN PTs ONLY in matplotlib
    grid_Info[10] = params[28]         # Plot LAGRANGIAN PTs + VELOCITY FIELD in matplotlib
    grid_Info[11] = params[29]         # Plot LAGRANGIAN PTs + VORTICITY colormap in matplotlib
    grid_Info[12] = params[30]         # Plot LAGRANGIAN PTs + MAGNITUDE OF VELOCITY 
                                       #     colormap in Matlab
    grid_Info[13] = params[31]         # Plot LAGRANGIAN PTs + PRESSURE colormap in matplotlib


    # MODEL STRUCTURE DATA STORED #
    model_Info = list(range(16))
    model_Info[0] = params[9]          # Springs: 0 (for no) or 1 (for yes) 
    model_Info[1] = params[10]         # Update_Springs: 0 (for no) or 1 (for yes)
    model_Info[2] = params[11]         # Target_Pts: 0 (for no) or 1 (for yes)
    model_Info[3] = params[12]         # Update_Target_Pts: 0 (for no) or 1 (for yes)
    model_Info[4] = params[13]         # Beams: 0 (for no) or 1 (for yes)
    model_Info[5] = params[14]         # Update_Beams: 0 (for no) or 1 (for yes)
    model_Info[6] = params[15]         # Muscle Activation (Length/Tension-Hill Model): 
                                       #     0 (for no) or 1 (for yes)
    model_Info[7] = params[16]         # Muscle Activation 3-ELEMENT HILL MODEL w/
                                       #     Length-Tension/Force-Velocity: 
                                       #     0 (for no) or 1 (for yes)
    model_Info[8] = params[17]         # Arbirtary External Force Onto Fluid Grid: 
                                       #     0 (for no) or 1 (for yes)
    model_Info[9] = params[18]         # Tracer Particles: 0 (for no) or 1 (for yes)
    model_Info[10]= params[19]         # Mass Points: 0 (for no) or 1 (for yes)
    model_Info[11]= params[20]         # Gravity: 0 (for no) or 1 (for yes)
    model_Info[12]= params[21]         # x-Component of Gravity vector
    model_Info[13]= params[22]         # y-Component of Gravity Vector
    model_Info[14]= params[23]         # Porous Media: 0 (for no) or 1 (for yes)
    model_Info[15]= params[24]         # Background Concentration Gradient: 
                                       #     0 (for no) or 1 (for yes)


    #-#-#-# DO THE IMMERSED BOUNDARY SOLVE!!!!!!!! #-#-#-#
    #[X, Y, U, V, xLags, yLags] = Driver.main(struct_name, mu, rho, grid_Info, dt, T_final, model_Info)
    
    #For debugging only!
    Driver.main(struct_name, mu, rho, grid_Info, dt, T_final, model_Info)
    
if __name__ == "__main__":
    main2d()