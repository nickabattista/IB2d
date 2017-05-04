'''-----------------------------------------------------------------------------------------------------------

 IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
 	fluid-structure interaction models. This version of the code is based off of
	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.

 Author: Nicholas A. Battista
 Email:  nick.battista@unc.edu
 Date Created: May 27th, 2015
 Institution: UNC-CH

 This code is capable of creating Lagrangian Structures using:
 	1. Springs
 	2. Beams (*torsional springs)
 	3. Target Points
   4. Mass Points
   5. Porous Points
	6. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-Tension)")
   7. 3-Element Hill Muscle Model

 One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting-lengths, etc
 
 There are a number of built in Examples, mostly used for teaching purposes. 
 
 If you would like us to add a specific muscle model, please let Nick (nick.battista@unc.edu) know.

-------------------------------------------------------------------------------------------------------------'''

import numpy as np
from numpy import pi as PI
import sys
import os

################################################################################
#
# FUNCTION: Removes comments from input file and creates temporary input file
#        
# Author:  Nick Battista, UNC-CH, nickabattista@gmail.com
#          Aaron Barrett, UNC-CH, abarret@live.unc.edu
#
################################################################################

def removeComments(input_file):

    FID = open(input_file)                 # Opens original input file
    input_file = input_file + '.temp'      # Makes temporary input file
    FID_temp = open(input_file, 'w')       # Creates/Opens for Writing temp. input file

    count = len(FID.readlines(  ))          # Counts # of lines in input file
    FID.seek(0)                             # Goes back to first line of file

    tline = FID.readline()
    num = 0                                 # Initialize counter
    while ( (isinstance(tline, str)) and (num < count+1) ):
        for i in range(len(tline)):         # Loops over all characters in line
            if( ('%' == tline[i]) or ( '#' == tline[i]) ):          # Checks if there is a comment and ends line there
                FID_temp.write('\n')        # MATLAB: fprintf(FID_temp,'\n')
                break                       # ends looping over that line
            
            FID_temp.write(tline[i])        # Writes character to temp input file 
        
        tline = FID.readline()              # Reads next line in
        num = num+1                         # Iterate counter

    FID.close                               # Closes original input file
    FID_temp.close                          # Closes temp input file

    return input_file


################################################################################
#
# FUNCTION: Stores the input file data as either a number or a string
#        
# Author:  Nick Battista, UNC-CH, nickabattista@gmail.com
#          Aaron Barrett, UNC-CH, abarret@live.unc.edu
#
################################################################################

def readValue(strInput):

    # FINDS INDICES OF QUOTATION MARKS, OR GIVES -1 FOR NUMBER
    index = 0
    while index < len(strInput):
        ind_aux = strInput.find('"', index)
        if ind_aux == -1:
            find_quote1 = -1
            break
        elif( index == 0 ):
            find_quote1 = ind_aux
            index = index + 1
        else:
            find_quote2 = index
            index = index+1

    if(find_quote1>-1):
        # Value is a string
        value = strInput[find_quote1+1:find_quote2]
        value.strip()
    else:
        # Value is a number
        value = float(strInput)

    return value


################################################################################
#
# FUNCTION: Reads in the input2d file information
#        
# Author:  Nick Battista, UNC-CH, nickabattista@gmail.com
#          Aaron Barrett, UNC-CH, abarret@live.unc.edu
#
################################################################################

def readBrace(FID):
    
    params = [[],[]]                    # Initialize params, MATLAB: params = cell(1,2)
    i = 1                               # Initialize index, i
    tline = FID.readline().rstrip() # Reads line in file FID, MATLAB: tline = fgets(fid)
    brace_close = tline.find('}')       # Finds index of character '} in line, 
                                                # MATLAB: brace_close = strfind(tline, '}')

    while( brace_close == -1 ):         # MATLAB: while(isempty(brace_close))
        
        brace_open = tline.find('{')                 # Finds character index of '{' in line
        if(brace_open>-1):                           # CHECK IF MORE SUB-CELLS ('lists' in python)
            params[0].append(tline[0:brace_open-1])  # MATLAB: params{i,1} = tline(1:brace_open(1)-1)
            f_inputs = readBrace(FID)
            params[1].extend(f_inputs)               # MATLAB: params{i,2} = readBrace(fid)
        
        else:
            find_equal = tline.find('=')             # Finds character index of '=' sign
            f_inputs = tline[0:find_equal-1].strip() # Reads all characters up to '=' sign and strips whitespace
            params[0].append(f_inputs)               # Stores all values up to equal sign in params[0]
                                                            # MATLAB: params{i,1} = tline(1:find_equal(1)-1)
            f_inputs = readValue(tline[find_equal+1:])
            

            params[1].append(f_inputs)               # Stores values after equal sign in params[1]
                                                            #MATLAB: params{i,2} = readValue(tline(find_equal(1)+1:end))

        tline = FID.readline().rstrip()                       # READ NEXT LINE
        brace_close = tline.find('}')                # Finds index of character '} in line, 
                                                            # MATLAB: brace_close = strfind(tline, '}')
        i = i+1                                      # Update iterator

    return params


################################################################################
#
# FUNCTION: Reads in the input2d file information if no braces (OLD FORMAT)
#        
# Author:  Nick Battista, UNC-CH, nickabattista@gmail.com
#          Aaron Barrett, UNC-CH, abarret@live.unc.edu
#
################################################################################

def readNoBrace(FID):

    params = [[],[]]                    # Initialize params, MATLAB: params = cell(1,2)
    num=0                               # Initialize counter, num
    tline = FID.readline()              # Reads line in file FID, MATLAB: tline = fgets(fid)
    count = len(FID.readlines())        # Counts # of lines in input file
    FID.seek(0)                         # Goes back to first line of file after counting

    while ( (isinstance(tline, str)) and (num<count+1) ): # MATLAB: while(ischar(tline))
        
        find_equal = tline.find('=')    # Finds character index of '=' sign        if(find_equal):

        params[0].append(tline[0:find_equal-1])  # Stores all values up to equal sign in params[0]
                                                            # MATLAB: params{i,1} = tline(1:find_equal(1)-1)
        f_inputs = readValue(tline[find_equal+1:])
        params[1].append(f_inputs)               # Stores values after equal sign in params[1]
                                                            # MATLAB: params{i,2} = readValue(tline(find_equal(1)+1:end))
        num=num+1;                               # Iterates num                                             

        tline = FID.readline()          # Reads next line in file FID, MATLAB: tline = fgets(fid)

    return params


################################################################################
#
# FUNCTION: TEST PRINTING OF READING IN FUNCTION
#        
# Author:  Nick Battista, UNC-CH, nickabattista@gmail.com
#
################################################################################
 
def please_Test_Reading(parameters):

    # PRINTS INPUT GROUP CATEGORY NAMES
    print(parameters[0][0][:]) # Prints 'Fluid Parameters'
    print(parameters[0][1][:]) # Prints 'Temporal Information'
    print(parameters[0][2][:]) # Prints 'Grid Parameters'
    print(parameters[0][:][:]) # Prints all category names

    # PRINTS PARAMETER NAME
    print(parameters[1][0][:]) # Prints ['mu',rho']
    print(parameters[1][0][0]) # Prints 'mu'
    print(parameters[1][0][1]) # Prints 'rho'

    # PRINTS PARAMETER VALUE
    print(parameters[1][1][:]) # Prints ['mu',rho'] = [0.01, 1.0] <-only prints numbers
    print(parameters[1][1][0]) # Prints 0.01
    print(parameters[1][1][1]) # Prints 1.0

    return 


################################################################################
#
# FUNCTION: Reads in the input2d file information
#        
# Author:  Nick Battista, UNC-CH, nickabattista@gmail.com
#          Aaron Barrett, UNC-CH, abarret@live.unc.edu
#
################################################################################

def please_Read_input2d_File(file_name):

    #file_name = 'input2d.IB2d'

    file_name = removeComments(file_name)    
    FID = open(file_name, 'r')              # Open file_name that's passed in
    parameters = [[],[]]                    # Create list ('cell-array')
    
    count = len(FID.readlines())            # Counts # of lines in input file
    FID.seek(0)                             # Goes back to first line of file after counting

    tline = FID.readline()                  # Read first line
    num = 0                                 # Initialize counter
    
    while ( (isinstance(tline, str)) and (num < count+1) ):

        find_brace = tline.find('{')        # Finds index for { in string, if not, gives -1
        if(find_brace>-1):
            parameters[0].append(tline[0:find_brace-1])  # MATLAB: parameters{i,1} = strtrim(tline(1:find_brace(1)-1))
            f_inputs = readBrace(FID)       # Returns list (cell array) with values between parenthesis
            parameters[1].append(f_inputs)  # Stores values after equal sign in params[1]
                                                         # MATLAB: parameters{i,2} = readBrace(FID)

        tline = FID.readline()              # Reads next line in
        num = num+1                         # Iterate counter

    # Prints out some information for testing purposes
    # please_Test_Reading(parameters)         

    FID.seek(0)                # Goes back to first line of file after counting

    # For reading in old input file format without braces (not necessary now)
    '''if (isempty(parameters{1,1})):
        # Whoops, no braces
        parameters = readNoBrace(FID)
    '''

    FID.close()                # Closes temp input file, MATLAB: close(FID)
    os.remove(file_name)       # Deletes temp input file, MATLAB: delete(file_name)

    return parameters          # Returns Parameters

#if __name__ == "__main__":
#    please_Read_input2d_File()




