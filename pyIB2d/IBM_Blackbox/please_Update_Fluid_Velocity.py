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
from numpy import pi as PI
from numpy import sin, cos
try:
    import pyfftw
    # create global variables
    fft2 = None
    ifft2 = None
    fft_mat = None
    ifft_mat = None
    FFTW = True
    print('Running with pyFFTW.')
except:
    print('Running without pyFFTW library.')
    FFTW = False
from Supp import D, DD

################################################################################
#
# This function solves the incompressible Navier-Stokes (NS) equations using 
#      Fast-Fourier Transform (FFT)
#      
#      x-Momentum Conservation: 
#           rho*u_t = -rho*u*u_x + rho*v*u_y + mu*laplacian(u) - p_x + Fx
#      y-Momentum Conservation: 
#           rho*v_t = -rho*u*v_x + rho*v*v_y + mu*laplacian(v) - p_y + Fy
#      Incompressibility:   u_x + v_y = 0.
#
################################################################################

def please_Update_Fluid_Velocity(U, V, Fx, Fy, rho, mu, grid_Info, dt, idX, idY):
    '''Fluid (Eulerian) Grid updated using Peskin's two-step algorithm, 
    where the advection terms are written in skew symmetric form.
    
    Args:
        U:         Eulerian grid of x-Velocities
        V:         Eulerian grid of y-Velocities
        Fx:        x-Forces arising from deformations in Lagrangian structure
        Fy:        y-Forces arising from deformatinos in Lagrangian structure
        rho:       Fluid density
        mu:        Fluid dynamic viscosity
        grid_Info: Vector of parameters relating to Eulerian grid, 
                       Lagrangian grid, numerical delta function
        dt:        Time-step
        idX/idY:   EULERIAN Index Matrices for FFT Operators
        
    Returns:
        U_h:
        V_h:
        U:
        V:
        p:'''

    # Initialize #
    Nx =   grid_Info[0]
    Ny =   grid_Info[1]
    Lx =   grid_Info[2]
    Ly =   grid_Info[3]
    dx =   grid_Info[4]
    dy =   grid_Info[5]
    supp = int(grid_Info[6])
    Nb =   grid_Info[7]
    ds =   grid_Info[8]


    # Construct EULERIAN Index Matrices
    #idX = np.tile(np.arange(Nx),(Nx,1))
    #idY = np.tile(np.arange(Ny),(Ny,1)).T
    
    # Create an FFTW plan if pyfftw is loaded
    if FFTW:
        global fft2, ifft2, fft_mat, ifft_mat
        if fft2 is None:
            fft_mat = pyfftw.empty_aligned((Ny,Nx), dtype='complex128')
            fft2 = pyfftw.builders.fft2(fft_mat)
            ifft_mat = pyfftw.empty_aligned((Ny,Nx), dtype='complex128')
            ifft2 = pyfftw.builders.ifft2(ifft_mat)
            

    # Create FFT Operator (used for both half time-step and full time-step computations)
    A_hat = 1 + 2*mu*dt/rho*( (sin(PI*idX/Nx)/dx)**2 + (sin(PI*idY/Ny)/dy)**2 )


    # # # # # # EVOLVE THE FLUID THROUGH HALF A TIME-STEP # # # # # # #


    # Compute 1ST and 2ND derivatives (using centered differencing)
    Ux = D(U,dx,'x')
    Uy = D(U,dy,'y')

    Uxx = DD(U,dx,'x')
    Uyy = DD(U,dy,'y')

    Vx = D(V,dx,'x')
    Vy = D(V,dy,'y')

    Vxx = DD(V,dx,'x')
    Vyy = DD(V,dy,'y')

    # Find derivatives of products U^2, V^2, and U*.
    U_sq = U**2
    V_sq = V**2
    UV = U*V

    U_sq_x = D(U_sq,dx,'x')
    V_sq_y = D(V_sq,dy,'y')

    UV_x = D(UV,dx,'x')
    UV_y = D(UV,dy,'y')

    #
    # Construct right hand side in linear system and perform FFT to
    #   take velocities to state space
    #
    rhs_u = give_RHS_HALF_Step_Velocity(dt,rho,Nx,Ny,U,Ux,Uy,U_sq_x,UV_y,V,Fx,'x')
    if FFTW:
        rhs_u_hat = np.array(fft2())
    rhs_v = give_RHS_HALF_Step_Velocity(dt,rho,Nx,Ny,V,Vx,Vy,V_sq_y,UV_x,U,Fy,'y')
    if FFTW:
        rhs_v_hat = fft2()
    if not FFTW:
        rhs_u_hat = np.fft.fft2(rhs_u)
        rhs_v_hat = np.fft.fft2(rhs_v)

    # Calculate Fluid Pressure
    p_hat = give_Fluid_Pressure(0.5*dt,rho,dx,dy,Nx,Ny,idX,idY,rhs_u_hat,rhs_v_hat)


    # Calculate Fluid Velocity and take IFFT to get velocities in real space
    u_hat = give_Me_Fluid_Velocity(0.5*dt,rho,dx,Nx,Ny,rhs_u_hat,p_hat,A_hat,idX,'x')
    if FFTW:
        U_h = np.real(np.array(ifft2(u_hat))) # the last here ifft2 replaces ifft_mat.
    v_hat = give_Me_Fluid_Velocity(0.5*dt,rho,dy,Nx,Ny,rhs_v_hat,p_hat,A_hat,idY,'y')
    if FFTW:
        V_h = np.real(np.array(ifft2()))
    if not FFTW:
        U_h = np.real(np.fft.ifft2(u_hat))   #Half-step velocity, u
        V_h = np.real(np.fft.ifft2(v_hat))   #Half-step velocity, v




    # # # # # # NOW EVOLVE THE FLUID THROUGH A FULL TIME-STEP # # # # # # #



    # Compute first derivatives (centred) at half step.
    U_h_x = D(U_h,dx,'x')
    U_h_y = D(U_h,dy,'y')
    V_h_x = D(V_h,dx,'x')
    V_h_y = D(V_h,dy,'y')

    # Computed derivatives of products U^2, V^2, and U*V at half step.
    U_h_sq = U_h**2
    V_h_sq = V_h**2
    U_h_V_h = U_h*V_h
    U_h_sq_x = D(U_h_sq,dx,'x')
    V_h_sq_y = D(V_h_sq,dy,'y')
    U_h_V_h_x = D(U_h_V_h,dx,'x')
    U_h_V_h_y = D(U_h_V_h,dy,'y')

    # Construct right hand side in linear system and peform FFT to take
    #   velocities to state space
    rhs_u = give_RHS_FULL_Step_Velocity(dt,mu,rho,Nx,Ny,U,U_h_x,U_h_y,
                                        U_h_sq_x,U_h_V_h_y,V,Fx,Uxx,Uyy,'x')
    if FFTW:
        rhs_u_hat = np.array(fft2()) #overwrites rhs_v_hat from before
    rhs_v = give_RHS_FULL_Step_Velocity(dt,mu,rho,Nx,Ny,V,V_h_x,V_h_y,
                                        V_h_sq_y,U_h_V_h_x,U,Fy,Vxx,Vyy,'y')
    if FFTW:
        rhs_v_hat = fft2()
    if not FFTW:
        rhs_u_hat = np.fft.fft2(rhs_u)
        rhs_v_hat = np.fft.fft2(rhs_v)  


    # Calculate Fluid Pressure
    p_hat  = give_Fluid_Pressure(dt,rho,dx,dy,Nx,Ny,idX,idY,rhs_u_hat,rhs_v_hat)
            

    # Calculate Fluid Velocity and take IFFT to get real velocities/pressure
    u_hat = give_Me_Fluid_Velocity(dt,rho,dx,Nx,Ny,rhs_u_hat,p_hat,A_hat,idX,'x')
    if FFTW:
        U = np.real(np.array(ifft2()))
    v_hat = give_Me_Fluid_Velocity(dt,rho,dy,Nx,Ny,rhs_v_hat,p_hat,A_hat,idY,'y')
    if FFTW:
        V = np.real(np.array(ifft2()))
        p = np.real(np.array(ifft2(p_hat)))
    if not FFTW:
        U = np.real(np.fft.ifft2(u_hat))
        V = np.real(np.fft.ifft2(v_hat))
        p = np.real(np.fft.ifft2(p_hat))

    return (U_h, V_h, U, V, p)



################################################################################
#
# FUNCTION: calculates the fluid velocity!
#
################################################################################

def give_Me_Fluid_Velocity(dt,rho,dj,Nx,Ny,rhs_VEL_hat,p_hat,A_hat,idMat,string):
    ''' Calculates the fluid velocity. Assigns result to FFTW plan if applicable.

    Args:
        dt: double
        rho: double
        dj: double
        Nx: int
        Ny: int
        rhs_VEL_hat: 2D array complex128
        p_hat: 2D array complex128
        A_hat: 2D array float
        idMat: 2D array int
        string: char
        
    Returns:
        vel_hat: fluid velocity'''

    if FFTW:
        global ifft_mat
        vel_hat = ifft_mat
    else:
        vel_hat = np.zeros((Ny,Nx),dtype='complex128') #initialize fluid velocity


    if string=='x':
        vel_hat[:,:] = ( rhs_VEL_hat - 1j*dt/(rho*dj)*sin(2*PI*idMat/Nx)*p_hat )/A_hat

    elif string=='y':
        vel_hat[:,:] = ( rhs_VEL_hat - 1j*dt/(rho*dj)*sin(2*PI*idMat/Ny)*p_hat )/A_hat
        
    return vel_hat

################################################################################
#
# FUNCTION: creates RHS with fluid velocity in FULL STEP computation
#
################################################################################

def give_RHS_FULL_Step_Velocity(dt,mu,rho,Nx,Ny,A,Ax,Ay,A_sq_j,AB_j,B,Fj,\
    Axx,Ayy,string):
    ''' Creates RHS with fluid velocity in FULL STEP computation.
        Assigns result to FFTW plan if applicable.
    
    Args:
        dt:
        mu:
        rho:
        Nx:
        Ny:
        A:
        Ax:
        Ay:
        A_sq_j:
        AB_j:
        B:
        Fj:
        Axx:
        Ayy:
        string:
        
    Returns:
        rhs:'''

    # Note: Fj -> j corresponds to either x or y.
    
    if FFTW:
        global fft_mat
        rhs = fft_mat #use matrix tied to fftw plan
    else:
        rhs = np.zeros((Ny,Nx)) #initialize rhs

    if string=='x':
        rhs[:,:] = A + dt/rho*( Fj + mu/2*(Axx+Ayy) - 0.5*rho*(A*Ax+B*Ay) -\
                .5*rho*(A_sq_j+AB_j) ) #RHS: u-component
        
    elif string=='y':
        rhs[:,:] = A + dt/rho*( Fj + mu/2*(Axx+Ayy) - 0.5*rho*(B*Ax+A*Ay) -\
                .5*rho*(AB_j+A_sq_j) ) #RHS: v-compoent
       
    return rhs




###########################################################################################
#
# FUNCTION: creates RHS with fluid velocity in HALF computation
#
###########################################################################################

def give_RHS_HALF_Step_Velocity(dt,rho,Nx,Ny,A,Ax,Ay,A_sq_j,AB_j,B,Fj,string):
    ''' Creates RHS with fluid velocity in HALF computation.
        Assigns result to FFTW plan if applicable.
    
    Args:
        dt:
        rho:
        Nx:
        Ny:
        A:
        Ax:
        Ay:
        A_sq_j:
        AB_j:
        B:
        Fj:
        string:
        
    Returns:
        rhs:'''

    # Note: Fj -> j corresponds to either x or y.

    if FFTW:
        global fft_mat
        rhs = fft_mat #use matrix tied to fftw plan
    else:
        rhs = np.zeros((Ny,Nx)) #initialize rhs

    if string=='x':
        #RHS: u-component
        rhs[:,:] = A + .5*dt/rho*( Fj - .5*rho*(A*Ax+B*Ay) - .5*rho*(A_sq_j+AB_j) )
        
    elif string=='y':
        #RHS: v-compoent
        rhs[:,:] = A + .5*dt/rho*( Fj - .5*rho*(B*Ax+A*Ay) - .5*rho*(AB_j+A_sq_j) )
       
    return rhs


################################################################################
#
# FUNCTION: calculates the fluid pressure
#
################################################################################

def give_Fluid_Pressure(dt,rho,dx,dy,Nx,Ny,idX,idY,rhs_u_hat,rhs_v_hat):
    ''' Calculates the fluid pressure. Result is FFTW friendly if applicable.
    
    Args:
        dt: double
        rho: double
        dx: double
        dy: double
        Nx: int
        Ny: int
        idX: 2D array of ints
        idY: 2D array of ints
        rhs_u_hat: 2D array of complex128
        rhs_v_hat: 2D array of complex128
        
    Returns:
        p_hat:'''

    if FFTW:
        p_hat = pyfftw.empty_aligned((Ny,Nx), dtype='complex128')
    else:
        p_hat = np.empty((Ny,Nx),dtype='complex128') #initialize fluid pressure
    
    num = -( 1j/dx*sin(2*PI*idX/Nx)*rhs_u_hat + 1j/dy*sin(2*PI*idY/Ny)*rhs_v_hat )
    den = ( dt/rho*( (sin(2*PI*idX/Nx)/dx)**2 + (sin(2*PI*idY/Ny)/dy)**2 ) )
    
    # Deal with nan in [0,0] entry
    num[0,0] = 0; den[0,0] = 1
    
    p_hat[:,:] = num/den

    # Zero out modes.
    p_hat[0,0] = 0
    p_hat[0,int(Nx/2)] = 0
    p_hat[int(Ny/2),int(Nx/2)] = 0
    p_hat[int(Ny/2),0] = 0

    return p_hat