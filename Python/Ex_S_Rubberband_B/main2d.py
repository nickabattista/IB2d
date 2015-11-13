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
    
    '''This is the "main" file, which ets called to run the Immersed Boundary Simulation. 
    It reads in all the parameters from "input2d", and sends them 
    off to the "IBM_Driver" function to actual perform the simulation.'''
    
    # READ-IN INPUT PARAMTERS #
    params,strName = give_Me_input2d_Parameters()
    
    