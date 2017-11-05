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
 
 For the Python port, I am going to throw a lot of supporting functions into
    here for convinence. That way they get loaded all at once, and are called
    by their name in an intuitive way. The functions (with their subfunctions
    following) are here in this order:
    -- please_Move_Lagrangian_Point_Positions
    -- give_NonZero_Delta_Indices_XY
    -- give_Eulerian_Lagrangian_Distance
    -- give_Delta_Kernel
    -- give_1D_NonZero_Delta_Indices
    -- please_Move_Massive_Boundary
    -- please_Update_Massive_Boundary_Velocity
    -- D
    -- DD
    -- please_Update_Adv_Diff_Concentration

----------------------------------------------------------------------------'''

import numpy as np
from math import sqrt
from numba import jit
    
################################################################################
#
# FUNCTION: Moves Lagrangian Point Positions by doing the integral,
#
#           " xLag_Next = xLag_Prev + dt* int( u(x,t) delta( x - xLag_n ) dX ) "
#
################################################################################

def please_Move_Lagrangian_Point_Positions(mu, u, v, xL_P, yL_P, xL_H, yL_H, x, y,\
    dt, grid_Info,porous_Yes,poroelastic_Yes,poroelastic_info,F_Poro):
    ''' Moves Lagrangian point positions
        u: 2D array
        v: 2D array
        xL_P:
        yL_P:
        xL_H:
        yL_H:
        x:
        y:
        dt:
        grid_Info:
        porous_Yes:
    
    Returns:
        xL_Next: 
        yL_Next:'''


    # Grid Info. grid_Info is a dict
    Nx =   grid_Info[0]
    Ny =   grid_Info[1]
    Lx =   grid_Info[2]
    Ly =   grid_Info[3]
    dx =   grid_Info[4]
    dy =   grid_Info[5]
    supp = int(grid_Info[6])
    Nb =   grid_Info[7]
    ds =   grid_Info[8]


    # Find indices where the delta-function kernels are non-zero for both x and y.
    xInds,yInds = give_NonZero_Delta_Indices_XY(xL_H, yL_H, Nx, Ny, dx, dy, supp)

    # ReSize the xL_H and yL_H matrices for use in the Dirac-delta function
    #        values to find distances between corresponding Eulerian data and them
    xLH_aux = xL_H % Lx
    yLH_aux = yL_H % Ly
    # Stack copies of the row vector and then transpose
    xL_H_ReSize = np.tile(xLH_aux,(supp**2,1)).T
    yL_H_ReSize = np.tile(yLH_aux,(supp**2,1)).T


    # Finds distance between specified Eulerian data and nearby Lagrangian data
    # x is a 1D array. x[xInds] is a 2D array of values in x
    if ( np.isscalar(xL_P) ):
        distX = give_Eulerian_Lagrangian_Distance(x[xInds], xL_H_ReSize, Lx)
        distY = give_Eulerian_Lagrangian_Distance(y[yInds], yL_H_ReSize, Ly)
    else:
        distX = give_Eulerian_Lagrangian_Distance(x[xInds], xL_H_ReSize, Lx)
        distY = give_Eulerian_Lagrangian_Distance(y[yInds], yL_H_ReSize, Ly)

    # Obtain the Dirac-delta function values.
    delta_X = give_Delta_Kernel(distX, dx)
    delta_Y = give_Delta_Kernel(distY, dy)

    # Perform Integral
    move_X, move_Y = give_Me_Perturbed_Distance(u,v,dx,dy,delta_X,delta_Y,xInds,yInds)

    # Update the Lagrangian Point Position.
    xL_Next = xL_P + (dt) * move_X
    yL_Next = yL_P + (dt) * move_Y

    # Update the Lagrangian Point Positions with poroelasticity.
    if poroelastic_Yes:
        #
        # poroelastic_info[:,1]: index of poroelastic point
        # poroelastic_info[:,2]: Brinkman constant
        #
        xL_Next[poroelastic_info[0:,0].astype(int)] = xL_Next[poroelastic_info[0:,0].astype(int)] + (  1/(mu*poroelastic_info[0:,1]) * F_Poro[0:,0] ) * dt
        yL_Next[poroelastic_info[0:,0].astype(int)] = yL_Next[poroelastic_info[0:,0].astype(int)] + (  1/(mu*poroelastic_info[0:,1]) * F_Poro[0:,1] ) * dt



    # Shift so that all values are in [0,Lx or Ly).
    if not porous_Yes:
        xL_Next = xL_Next % Lx
        yL_Next = yL_Next % Ly

    return (xL_Next, yL_Next)



################################################################################
# FUNCTION: Computes the integral to move each Lagrangian Pt!
################################################################################

def give_Me_Perturbed_Distance(u,v,dx,dy,delta_X,delta_Y,xInds,yInds):
    ''' Computes the integral to move each Lagrangian Pt.
    
    Args:
        u:        x-component of velocity (2D array)
        v:        y-component of velocity (2D array)
        delta_X:  values of Dirac-delta function in x-direction (2D array)
        delta_Y:  values of Dirac-delta function in y-direction (2D array)
        xInds:    x-Indices on fluid grid (2D array)
        yInds:    y-Indices on fluid grid (2D array)'''

    # Compute integrand 'stencil' of velocity x delta for each Lagrangian Pt!
    # Fancy indexing allows us to do this directly
    mat_X = u[yInds,xInds]*delta_X*delta_Y
    mat_Y = v[yInds,xInds]*delta_X*delta_Y

            
    # Approximate Integral of Velocity x Delta for each Lagrangian Pt!
    move_X = mat_X.sum(1) * (dx*dy)
    move_Y = mat_Y.sum(1) * (dx*dy)

    return (move_X, move_Y)




############################################################################################
#
# FUNCTION: finds the indices on the Eulerian grid where the 1D Dirac-delta
# kernel is possibly non-zero in BOTH (x,y) directions
#
############################################################################################

def give_NonZero_Delta_Indices_XY(xLag, yLag, Nx, Ny, dx, dy, supp):
    ''' Find indices where 1D Dirac-delta kernel is non-zero in both (x,y)
    
    Args:
        xLag: gives x-coordinate of Lagrangian position
        yLag: gives y-coordinate of Lagrangian position
        Nx:   # of Eulerian grid pts. in x-dimension
        Ny:   # of Eulerian grid pts. in y-dimension
        dx:   spatial-step along x-dimension of Eulerian grid
        dy:   spatial-step along y-dimension of Eulerian grid
        supp: size of support of the Dirac-delta kernel (should be even)
        
    Returns:
        xInds: x index
        yInds: y index'''


    #Give x-dimension Non-Zero Delta Indices
    xIndsAux = give_1D_NonZero_Delta_Indices(xLag, Nx, dx, supp)

    #Repeat x-Indices for Non-Zero y-Indices!
    xInds = []
    #Sets up x-INDEX matrix bc we consider BOTH dimensions
    xInds = np.tile(xIndsAux,(1,supp)) #tiles matrix in horiz direction
    

    #Give y-dimension Non-Zero Delta Indices
    yIndsAux = give_1D_NonZero_Delta_Indices(yLag, Ny, dy, supp)

    #Repeat y-Indices for Non-Zero x-Indices!
    yInds = np.repeat(yIndsAux,supp,axis=1) #repeats each element horizontally
                                            #   supp number of times
    #Sets up y-INDEX matrix bc we consider BOTH dimensions
    
    #these are indices, so return ints
    return (xInds.astype('int'),yInds.astype('int'))



################################################################################
#
# FUNCTION distance between Eulerian grid data, x, and Lagrangian grid data, y, 
#          at specifed pts typically and makes sure the distances are 
#          correct for a periodic [0,L] domain.
#
################################################################################

def give_Eulerian_Lagrangian_Distance(x, y, L):
    ''' Find dist. between Eulerian grid data and Lagrangian grid data.
    [0,L] has periodic boundary condition, so in actuality, the greatest
    distance possible is L/2.
    
    Args:
        x,y: two matrices that you find the distance between
            (x-typically Eulerian data, y-typically Lagrangian data)
        L: length of domain, i.e., [0,L]
        
    Returns:
        distance: distance'''

    distance = abs( x - y )
    distance = np.minimum(distance,L-distance) #element-wise minima
    
    return distance



###########################################################################
#
# FUNCTION: computes a discrete approx. to a 1D Dirac-delta function over a
# specified matrix, x, and spatial step-size, dx. It will have support in
# [x-2dx, x+2dx]
#
###########################################################################

@jit(nopython=True)
def give_Delta_Kernel(x,dx):
    ''' Computes discrete approx. to 1D delta func over x in [x-2dx,x+2dx].
    
    Args:
        x:  Values in which the delta function will be evaulated
        dx: Spatial step-size of grid
        
    Returns:
        delta: delta function with support [x-2dx,x+2dx]'''

    # Computes Dirac-delta Approximation.
    RMAT = np.abs(x)/dx

    #Alias the data for cleaner writing of the following step
    #   RMAT is altered, but it will not be reused.
    delta = RMAT

    #Loops over to find delta approximation
    row,col = x.shape
    for ii in range(row):
        for jj in range(col):
            
            r = RMAT[ii,jj]
            
            # Approximate Discrete Delta Function
            #if r <= 2:
            #    delta[ii,jj] = 0.25*( 1 + cos(pi*r/2) )    
            #else:
            #    delta[ii,jj] = 0        

            # PESKIN 4-Point Discrete Delta Function
            if r<1:
                delta[ii,jj] = ( (3 - 2*r + sqrt(1 + 4*r - 4*r*r) ) / (8*dx) )
            elif (r<2) and (r>=1):
                delta[ii,jj] = ( (5 - 2*r - sqrt(-7 + 12*r - 4*r*r) ) / (8*dx) )
            else:
                delta[ii,jj] = 0
    return delta



###########################################################################
#
# FUNCTION finds the indices on the Eulerian grid where the 1D Dirac-delta
# kernel is possibly non-zero is x-dimension.
#
###########################################################################

def give_1D_NonZero_Delta_Indices(lagPts_j, N, dx, supp):
    ''' Find the indices on Eulerian grid where 1D delta is non-zero in x dim.
    
    Args:
        lagPts_j: row of lagrangian pts for specific coordinate, j= x or y.
        N:        # spatial resolution of Eulerian grid in each dimension
        dx:       Spatial step-size on Eulerian (fluid) grid
        supp:     Size of support of the Dirac-delta kernel (should be even)
        
    Returns:
        indices'''


    # Finds the index of the lower left Eulerian pt. to Lagrangian pt..
    ind_Aux = np.floor(lagPts_j/dx + 1)

    # Get all the different x indices that must be considered.
    # ind_Aux is 1D. Create 2D array with supp # of columns of ind_Aux
    indices = np.tile(ind_Aux,(supp,1)).T #stack row vectors then transpose

    #
    indices += -supp/2+1+np.arange(supp) #arange returns row array, which
                                         # broadcasts down each column.

    # Translate indices between {0,2,..,N-1}
    indices = (indices-1) % N
    
    return indices
    
################################################################################
#
# FUNCTION: update 'massive' immersed boundary position
#
################################################################################

def please_Move_Massive_Boundary(dt_step,mass_info,mVelocity):
    ''' Update 'massive' immersed boundary position
    
    Args:
        dt_step: desired time-step for this position
        mass_info: col 1: lag index for mass pt
                   col 2: massive x-Lag Value
                   col 3: massive y-Lag Value
                   col 4: 'mass-spring' stiffness parameter
                   col 5: MASS parameter value
        mVelocity  col 1: x-directed Lagrangian velocity
                   col 2: y-directed Lagrangian velocity

    Returns:
        mass_info:
        massLagsOld:'''

    massLagsOld = np.array(mass_info[:,(1, 2)])

    # update x-Positions
    mass_info[:,1] = mass_info[:,1] + dt_step*mVelocity[:,0]

    # update y-Positions
    mass_info[:,2] = mass_info[:,2] + dt_step*mVelocity[:,1]
    
    return (mass_info, massLagsOld)
    
    
    
############################################################################################
#
# FUNCTION: update 'massive' immersed boundary velocity
#
############################################################################################

def please_Update_Massive_Boundary_Velocity(dt_step,mass_info,mVelocity,\
    F_Mass_Bnd,gravity_Info):
    ''' Update 'massive' immersed boundary velocity
    
    Args:
        dt_step: desired time-step for this position
        mass_info:   col 1: lag index for mass pt
                     col 2: massive x-Lag Value
                     col 3: massive y-Lag Value
                     col 4: 'mass-spring' stiffness parameter
                     col 5: MASS parameter value
        mVelocity    col 1: x-directed Lagrangian velocity
                     col 2: y-directed Lagrangian velocity
        F_Mass_Bnd   col 1: x-directed Lagrangian force
                     col 2: y-directed Lagrangian force
        gravity_Info col 1: flag if considering gravity (0 = NO, 1 = YES)
                     col 2: x-component of gravity vector (NORMALIZED PREVIOUSLY)
                     col 3: y-component of gravity vector (NORMALIZED PREVIOUSLY)
                     
    Returns:
        mVelocity_h:'''

    ids = mass_info[:,0].astype('int')
    mVelocity_h = np.empty(mVelocity.shape)

    if gravity_Info[0] == 1:
         
        g = 9.80665 #m/s^2
        
        # update x-Velocity
        mVelocity_h[:,0] = mVelocity[:,0] - dt_step * \
        ( F_Mass_Bnd[ids,0]/mass_info[:,4] - g*gravity_Info[1] )

        # update y-Velocity
        mVelocity_h[:,1] = mVelocity[:,1] - dt_step * \
        ( F_Mass_Bnd[ids,1]/mass_info[:,4] - g*gravity_Info[2] )
        
    else:
        
        # update x-Velocity
        mVelocity_h[:,0] = mVelocity[:,0] - dt_step*F_Mass_Bnd[ids,0]/mass_info[:,4]

        # update y-Velocity
        mVelocity_h[:,1] = mVelocity[:,1] - dt_step*F_Mass_Bnd[ids,1]/mass_info[:,4]
        
    return mVelocity_h



########################################################################
#
# FUNCTION: Finds CENTERED finite difference approximation to 1ST
# Derivative in specified direction by input, dz, and 'string'. 
# Note: It automatically accounts for periodicity of the domain.
#
########################################################################

def D(u,dz,string):
    ''' Finds centered 1st derivative in specified direction
    
    Args:
        u:      velocity 
        dz:     spatial step in "z"-direction
        string: specifies which 1ST derivative to take (to enforce periodicity)
        
    Returns:
        u_z:'''

    
    u_z = np.zeros(u.shape)

    if string=='x':

        length = u.shape[1]      # number of pts along X direction

        #For periodicity on ends
        u_z[:,0] = ( u[:,1] - u[:,length-1] ) / (2*dz)
        u_z[:,-1]= ( u[:,0] - u[:,length-2] ) / (2*dz)

        #Standard Centered Difference
        u_z[:,1:length-1] = ( u[:,2:length] - u[:,:length-2] ) / (2*dz)

    elif string=='y':
        
        length = u.shape[0]      # number of pts along Y direction

        #For periodicity on ends
        u_z[0,:] = ( u[1,:] - u[length-1,:] ) / (2*dz)
        u_z[length-1,:] = ( u[0,:] - u[length-2,:] ) / (2*dz)

        #Standard Centered Difference
        u_z[1:length-1,:] = ( u[2:length,:] - u[:length-2,:] ) / (2*dz)
        
    else:
        
        print('\n\n\n ERROR IN FUNCTION D FOR COMPUTING 1ST DERIVATIVE\n')
        print('Need to specify which desired derivative, x or y.\n\n\n')
           
    return u_z


########################################################################
#
# FUNCTION: Finds CENTERED finite difference approximation to 2ND
# DERIVATIVE in z direction, specified by input and 'string' 
# Note: It automatically accounts for periodicity of the domain.
#
########################################################################

def DD(u,dz,string):
    ''' Finds centered 2nd derivative in z direction, specified by input & string
    
    Args:
        u:      velocity 
        dz:     spatial step in "z"-direction
        string: specifies which 2ND derivative to take (to enforce periodicity)
        
    Returns:
        u_zz:'''

    u_zz = np.zeros(u.shape)

    if string=='x':

        length = u.shape[1]      # number of pts along X direction


        #For periodicity on ends
        u_zz[:,0] =  ( u[:,1] - 2*u[:,0]   + u[:,length-1] )   / (dz**2)
        u_zz[:,-1] = ( u[:,0] - 2*u[:,length-1] + u[:,length-2] ) / (dz**2)

        #Standard Upwind Scheme (Centered Difference)
        u_zz[:,1:length-1] = (u[:,2:length]-2*u[:,1:length-1]+u[:,:length-2])\
                                / (dz**2)

    elif string=='y':

        length = u.shape[0]      # number of pts along Y direction

        #For periodicity on ends
        u_zz[0,:] =  ( u[1,:] - 2*u[0,:]   + u[length-1,:] )   / (dz**2)
        u_zz[-1,:]= ( u[0,:] - 2*u[length-1,:] + u[length-2,:] ) / (dz**2)

        #Standard Upwind Scheme (Centered Difference)
        u_zz[1:length-1,:] = (u[2:length,:]-2*u[1:length-1,:]+u[:length-2,:])\
                                / (dz**2)

    else:
        
        print('\n\n\n ERROR IN FUNCTION DD FOR COMPUTING 2ND DERIVATIVE\n')
        print('Need to specify which desired derivative, x or y.\n\n\n')
        
    return u_zz



###########################################################################
#
# def: Setting up advection-diffusion solver
#
###########################################################################

def please_Update_Adv_Diff_Concentration(C,dt,dx,dy,uX,uY,k):

    # C:     concentration 
    # dt:    time-step
    # dx,dy: spatial steps in x and y, respectively
    # uX:    x-Component of Velocity
    # uY:    y-Component of Velocity
    # k:     diffusion coefficient


    # Performs Upwind Advection WITHOUT Time-Splitting
    #C = perform_Time_noSplit_Upwind(C,dt,dx,dy,uX,uY,k)

    # Performs Upwind Advection w/ Time-Splitting
    C = perform_Time_Split_Upwind(C,dt,dx,dy,uX,uY,k)

    #laplacian_C=1 # DUMMY VARIABLE (laplacian not used anywhere else in code.)

    return C

###########################################################################
#
# def: Advection-Diffusion Split Upwind Method
#
###########################################################################

def perform_Time_noSplit_Upwind(C,dt,dx,dy,uX,uY,k):

    # Compute Necessary Derivatives (Note: these calculations could be parallalized)
    Cx = give_Necessary_Derivative(C,dx,uX,'x')
    Cy = give_Necessary_Derivative(C,dy,uY,'y') 
    Cxx = DD(C,dx,'x')
    Cyy = DD(C,dy,'y')

    # Forms Laplacian
    laplacian_C = Cxx+Cyy

    # UPWIND
    C = C + dt * ( k*(laplacian_C) - uX*Cx - uY*Cy )

    return C

###########################################################################
#
# def: Advection-Diffusion Split Upwind Method
#
###########################################################################

def perform_Time_Split_Upwind(C,dt,dx,dy,uX,uY,k):

    # Compute Necessary Derivatives for x-Advection 
    Cx = give_Necessary_Derivative(C,dx,uX,'x')
    Cxx = DD(C,dx,'x')

    # Time-step #1 (give auxillary)
    C = C + dt * ( k*(Cxx) - uX*Cx )

    # Compute Necessary Derivatives for y-Advection 
    Cy = give_Necessary_Derivative(C,dy,uY,'y') 
    Cyy = DD(C,dy,'y')

    # Time-step #2 (give next iteration)
    C = C + dt * ( k*(Cyy) - uY*Cy )

    return C


###########################################################################
#
# def: Computes derivative based on sign of Velocity, u, using UPWIND
# approach
#
###########################################################################

def give_Necessary_Derivative(C,dz,uZ,string):

    C_z = np.zeros(C.shape)
    signs = np.sign(uZ)
    [lenY,lenX] = uZ.shape

    if string=='x':

        #For periodicity on ends w/ UPWIND
        for i in range(0,lenY):

            #left side of grid
            if signs[i,0] <= 0: 
                C_z[i,0] =  ( C[i,1] - C[i,0] ) / (dz)
            else:
                C_z[i,0] =  ( C[i,0] - C[i,lenX-1] ) / (dz)

            #right side of grid
            if signs[i,lenX-1] <= 0: 
                C_z[i,lenX-1] =  ( C[i,0] - C[i,lenX-1] ) / (dz)
            else:
                C_z[i,lenX-1] =  ( C[i,lenX-1] - C[i,lenX-2] ) / (dz)


        #Standard Upwind 
        for i in range(0,lenY):
            for j in range(1,lenX-2):
                if signs[i,j] <= 0:
                    C_z[i,j] = ( C[i,j+1] - C[i,j] ) / (dz)
                else:
                    C_z[i,j] = ( C[i,j] - C[i,j-1] ) / (dz)

        # Ends x-Direction calculation #

    elif string=='y':

        #For periodicity on ends w/ UPWIND
        for i in range(0,lenX):

            #bottom of grid
            if signs[0,i] <= 0: 
                C_z[0,i] =  ( C[1,i] - C[0,i] ) / (dz)
            else:
                C_z[0,i] =  ( C[0,i] - C[lenY-1,i] ) / (dz)

            #top of grid
            if signs[lenY-1,i] <= 0: 
                C_z[lenY-1,i] =  ( C[0,i] - C[lenY-1,i] ) / (dz)
            else:
                C_z[lenY-1,i] =  ( C[lenY-1,i] - C[lenY-2,i] ) / (dz)

        #Standard Upwind
        for i in range(1,lenY-2):
            for j in range(0,lenX):
                if signs[i,j] <= 0:
                    C_z[i,j] = ( C[i+1,j] - C[i,j] ) / (dz)
                else:
                    C_z[i,j] = ( C[i,j] - C[i-1,j] ) / (dz)


        # Ends y-Direction calculation #

    else:

        print('\n\n\n ERROR IN def D FOR COMPUTING 1ST DERIVATIVE\n')
        print('Need to specify which desired derivative, x or y.\n\n\n') 

    
    return C_z



###########################################################################
#
# FUNCTION: Setting up advection-diffusion solver
#
###########################################################################

#def please_Update_Adv_Diff_Concentration(C,dt,dx,dy,uX,uY,k):
    '''Setting up advection-diffusion solver
    
    Note: This function alters C internally!
    
    Args:
        C:     concentration 
        dt:    time-step
        dx,dy: spatial steps in x and y, respectively
        uX:    x-Component of Velocity
        uY:    y-Component of Velocity
        k:     diffusion coefficient
        
    Returns:
        C:'''

  #  # Compute Necessary Derivatives for x-Advection 
  #  Cx = D(C,dx,'x')
  #  Cxx = DD(C,dx,'x')
   
  #  # Time-step #1 (give auxillary) 
  #  C = C + dt * ( k*(Cxx) - uX*Cx )

  #  # Compute Necessary Derivatives for y-Advection 
  #  Cy = D(C,dy,'y')
  #  Cyy = DD(C,dy,'y')

  #  # Time-step #2 (give next iterative for C) 
  #  C = C + dt * ( k*(Cyy) - uY*Cy )

  #  return C

