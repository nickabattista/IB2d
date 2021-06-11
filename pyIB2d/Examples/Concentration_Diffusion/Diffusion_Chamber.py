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
from math import cos, sin, pi, sqrt
import matplotlib.pyplot as plt

################################################################################
#
# FUNCTION: creates the CHANNEL_CHANNEL-EXAMPLE geometry 
#    and prints associated input files
#
################################################################################

def Diffusion_Chamber():

    #
    # Grid Parameters (MAKE SURE MATCHES IN input2d !!!)
    #
    Nx =  64        # # of Eulerian Grid Pts. in x-Direction (MUST BE EVEN!!!)
    Ny =  64        # # of Eulerian Grid Pts. in y-Direction (MUST BE EVEN!!!)
    Lx = 1.0        # Length of Eulerian Grid in x-Direction
    Ly = 1.0        # Length of Eulerian Grid in y-Direction
    dx = Lx/Nx      # Spatial step size in x-Direction
    dy = Ly/Ny      # Spatial step size in y-Direction


    # Immersed Structure Geometric / Dynamic Parameters #
    ds= min(Lx/(2*Nx),Ly/(2*Ny))  # Lagrangian spacing
    struct_name = 'chamber' # Name for .vertex, .spring, etc files.


    # Call function to construct geometry
    xLag,yLag = give_Me_Immsersed_Boundary_Geometry(ds,Lx,Ly)


    # Plot Geometry to test
    plt.plot(xLag,yLag,'r-')
    plt.plot(xLag,yLag,'*')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show(block=True)


    #
    # GIVES INITIAL CONCENTRATION 
    #
    C = give_Me_Initial_Concentration(Lx,Ly,Nx,Ny,dx,dy)



    # Prints .vertex file!
    print_Lagrangian_Vertices(xLag,yLag,struct_name)


    # Prints .concentration file!
    kDiffusion = 0.025
    print_Concentration_Info(Nx,Ny,C,kDiffusion,struct_name)


    # Prints .spring file!
    #k_Spring = 1e7
    #print_Lagrangian_Springs(xLag,yLag,k_Spring,ds_Rest,struct_name)


    # Prints .beam file!
    #k_Beam = 0.5; C = 0.0
    #print_Lagrangian_Beams(xLag,yLag,k_Beam,C,struct_name)


    # Prints .target file!
    k_Target = 1e7
    print_Lagrangian_Target_Pts(xLag,k_Target,struct_name)

########################################################################
#
# FUNCTION: prints VERTEX points to a file called rubberband.vertex
#
########################################################################

def print_Lagrangian_Vertices(xLag,yLag,struct_name):

    N = len(xLag)
    
    with open(struct_name + '.vertex', 'w') as vertex_fid:
        vertex_fid.write('{0}\n'.format(N))

        #Loops over all Lagrangian Pts.
        for s in range(N):
            X_v = xLag[s]
            Y_v = yLag[s]
            vertex_fid.write('{0:1.16e} {1:1.16e}\n'.format(X_v, Y_v))


########################################################################
#
# FUNCTION: prints CONCENTRATION INFO to file called
#           'struct_name'.concentration
#
########################################################################
    
def print_Concentration_Info(Nx,Ny,C,kDiffusion,struct_name):

    with open(struct_name + '.concentration','w') as con_fid:
        con_fid.write('{0}\n'.format(kDiffusion))
        
        for ii in range(Ny):
            for jj in range(Nx):
                con_fid.write('{0:1.16e} '.format(C[ii,jj]))
            con_fid.write('\n') 
    
    
########################################################################
#
# FUNCTION: prints Vertex points to a file called rubberband.vertex
#
########################################################################

def print_Lagrangian_Target_Pts(xLag,k_Target,struct_name):

    N = len(xLag)
    
    with open(struct_name + '.target','w') as target_fid:
        target_fid.write('{0}\n'.format(N))
        
        #Loops over all Lagrangian Pts.
        for s in range(N):
            target_fid.write('{0} {1:1.16e}\n'.format(s,k_Target))

    
    
################################################################################
#
# FUNCTION: prints BEAM (Torsional Spring) points to a file called rubberband.beam
#
################################################################################

def print_Lagrangian_Beams(xLag,yLag,k_Beam,C,struct_name):

    # k_Beam: beam stiffness
    # C: beam curvature
    
    N = len(xLag) # NOTE: Total number of beams = Number of Total Lag Pts. - 2
    
    with open(struct_name + '.beam','w') as beam_fid:
        beam_fid.write('{0}\n'.format(N))

        #spring_force = kappa_spring*ds/(ds**2)

        #BEAMS BETWEEN VERTICES
        for s in range(1,N-1):
            if s <= N-2:
                beam_fid.write('{0} {1} {2} {3:1.16e} {4:1.16e}\n'.format(
                    s-1,s,s+1,k_Beam,C))
            else:
                #Case s==N-1 (this does not occur)
                beam_fid.write('{0} {1} {2} {3:1.16e} {4:1.16e}\n'.format(
                    s-1,s,0,k_Beam,C))


########################################################################
#
# FUNCTION: prints SPRING points to a file called rubberband.spring
#
########################################################################

def print_Lagrangian_Springs(xLag,yLag,k_Spring,ds_Rest,struct_name):

    N = len(xLag)
    
    with open(struct_name + '.spring','w') as spring_fid:
        spring_fid.write('{0}\n'.format(N))
        
        #spring_force = kappa_spring*ds/(ds**2)

        #SPRINGS BETWEEN VERTICES
        for s in range(N):
            if s < N-1:
                spring_fid.write('{0} {1} {2:1.16e} {3:1.16e}\n'.format(
                    s,s+1,k_Spring,ds_Rest))
            else:
                #Case s==N-1
                spring_fid.write('{0} {1} {2:1.16e} {3:1.16e}\n'.format(
                    s,0,k_Spring,ds_Rest))
    
    

########################################################################
#
# FUNCTION: creates the Lagrangian structure geometry
#
########################################################################

def give_Me_Immsersed_Boundary_Geometry(ds,Lx,Ly):
    
    # The immsersed structure is a channel #
    xH = list(np.arange(2*ds,Lx-2*ds+ds/100,ds))  #xPts on horizontal part
    yHB = 2*ds*np.ones(len(xH))               #yPts on BOTTOM horizontal part
    yHT = Ly - yHB                            #yPts on TOP horizontal part
    
    yV = list(np.arange(3*ds,Ly-3*ds+ds/100,ds))  #xPts on horizontal part
    xVL = 2*ds*np.ones(len(yV))               #xPts on LEFT vertical part
    xVR = Lx - xVL                            #xPts on RIGHT vertical part
    
    #ORDER COMPONENTS CCW
    xLag = np.concatenate((xH,xVR,xH[::-1],xVL))
    yLag = np.concatenate((yHB,yV,yHT,yV[::-1]))
    
    return (xLag,yLag)

   
###################################################################
#
# FUNCTION: gives initial concentration gradient inside channel
#
###################################################################

def give_Me_Initial_Concentration(Lx,Ly,Nx,Ny,dx,dy):
    
    #WHERE OUTER TUBE LIES
    xMin = 0.475
    xMax = 0.525
    yMin = 0.475
    yMax = 0.525
    x = list(np.arange(0,Lx+dx/100,dx))
    y = list(np.arange(0,Ly+dy/100,dy))
    inds = give_Me_Indices_To_Apply_Force(x,y,xMin,xMax,yMin,yMax)
    
    C = np.zeros((Ny,Nx))
    for ii in range(inds.shape[0]):
        C[inds[ii,0],inds[ii,1]] = 1.0

    return C


######################################################################
#
# FUNCTION: computes indices for placing initial concentration 
#
#######################################################################

def give_Me_Indices_To_Apply_Force(x,y,xMin,xMax,yMin,yMax):
    '''There's almost certainly a better way to do this... argmax?'''
    
    
    jj=0
    noMinYet = True
    while noMinYet:
        if x[jj] >= xMin:
            iX_min = jj
            noMinYet = False
        jj += 1
        
    jj=len(x)-1
    noMaxYet = True
    while noMaxYet:
        if x[jj] <= xMax:
            iX_max = jj
            noMaxYet = False
        jj -= 1
        
    jj=0
    noMinYet = True
    while noMinYet:
        if y[jj] >= yMin:
            iY_min = jj
            noMinYet = False
        jj += 1

    jj=len(y)-1
    noMaxYet = True
    while noMaxYet:
        if y[jj] <= yMax:
            iY_max = jj
            noMaxYet = False
        jj -= 1
        
    iX_Vec = list(range(iX_min,iX_max+1))
    iY_Vec = list(range(iY_min,iY_max+1))
    
    n=0
    inds = np.zeros((len(iX_Vec)*len(iY_Vec),2))
    for ii in range(len(iX_Vec)):
        for jj in range(len(iY_Vec)):
            inds[n,0] = iX_Vec[ii]
            inds[n,1] = iY_Vec[jj]
            n += 1
            
    return inds
    
if __name__ == "__main__":
    Diffusion_Chamber()