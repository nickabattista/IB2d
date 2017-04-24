from please_Read_input2d_File import *

def test():

    # Read in parameters
    params = please_Read_input2d_File('input2d.IB2d')

    please_Test_Reading(params)

     # PRINTS INPUT GROUP CATEGORY NAMES
    #print(parameters[0][0][:]) # Prints 'Fluid Parameters'
    #print(parameters[0][1][:]) # Prints 'Temporal Information'
    #print(parameters[0][2][:]) # Prints 'Grid Parameters'
    #print(parameters[0][:][:]) # Prints all category names

    # PRINTS PARAMETER NAME
    #print(parameters[1][0][:]) # Prints ['mu',rho']
    #print(parameters[1][0][0]) # Prints 'mu'
    #print(parameters[1][0][1]) # Prints 'rho'

    # PRINTS PARAMETER VALUE
    #print(parameters[1][1][:]) # Prints ['mu',rho'] = [0.01, 1.0] <-only prints numbers
    #print(parameters[1][1][0]) # Prints 0.01
    #print(parameters[1][1][1]) # Prints 1.0
   
    # params[i][j]:  i -> group name OR 2x2 lists
    #                   j -> either string parameter name or variable
    #                   k -> gives 

    #print(params)
    #print(params[0][:][:].index('Grid_Parameters'))
    ind = params[0][:].index('Grid_Parameters')
    grid_Info = params[1][ind][:]
    
    print(grid_Info)
    print(grid_Info[0][0]) # Prints Nx
    print(grid_Info[0][1]) # Prints Ny
    print(grid_Info[1][0])
    print(grid_Info[1][1])

    ind2 = grid_Info[0][:].index('Nx')
    print(grid_Info[0][ind2])
    print(grid_Info[1][ind2])

    #t = grid_Info2[1]+4
    #print(grid_Info)
    #print(grid_Info2)
    #print(grid_Info2[1])
    #print(t)

if __name__ == "__main__":
    test()
    