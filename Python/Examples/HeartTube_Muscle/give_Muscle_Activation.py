##############################################################################
#
# FUNCTION: Returns vector of muscle activation forces for a single Lag.Pt.
#
##############################################################################

import numpy as np
from math import cos, sin, pi, sqrt, exp

def give_Muscle_Activation(v,LF,LFO,SK,a,b,Fmax,current_time,xPt,xLag):

    # current_time: current time in simulation (s)
    # xLag: vector of all x-Lagrangian Pts

    # Fm = a_f * Fmax * F1(Lf) * F2(Vf) 
    #    = a_f * Fmax *exp( -( (Q-1)/SK )^2 ) * (1/P0)*(b*P0-a*v)/(v+b); Q = LF/LFO
    #

    # a_f: coefficient that triggers contraction (traveling square or Gaussian wave?)

    # Length Tension Model Parameters (F1)
    # Fmax: maximum isometric force produced at the optimum length of the muscle fibers
    # LF:   length of the muscle fibers
    # LFO:  length at which the muscle fibers exert their maximum tension
    # SK:   constant specific for each muscle where SK > 0.28.


    # Hill Model Parameters (F2)
    # P0:   maximum load w/ NO contraction
    # a:    
    # b:    
    # v:    velocity of muscle expansion/contraction


    # Length Tension Model Parameters #
    Q = LF/LFO
    F1 = np.exp(-((Q-1)/SK)**2)

    # Hill Model #
    P0 = Fmax   #Same as Fmax
    F2 = (1/P0)*(b*P0-a*v)/(v+b)

    # Get Activation Coefficient #
    af_Val = give_Traveling_Triggering_Coefficient(current_time,xLag,xPt)

    # Actually Compute Muscle Force #
    Fm = af_Val*Fmax*F1*F2
    
    return Fm

#########################################################################
#
# FUNCTION: Returns value of Activation Trigger at specific time
#
#########################################################################


def give_Traveling_Triggering_Coefficient(current_time,xLag,xPt):

    # current_time: current_time in simulation
    # xLag:         x-Lagrangian pts (both bottom and top of tube)
    # xPt: x-Pt of interest

    t = current_time                # current time
    freq = 10                       # frequency of traveling wave down tube
    Ltube = xLag[len(xLag)//2 - 1] - xLag[0]  # length of heart tube
    buff = 0.25                     # percent-buff on each end of heart tube
    L_AR = Ltube*(1-2*buff)         # length of the actual activation region of tube
    SqWidth = L_AR/10               # width of the square traveling activation wave 
    v = (L_AR-SqWidth) * freq       # traveling wave velocity
    k = 2*pi*freq/v                 # Wave-number for traveling wave

    t = np.fmod(t,1/freq)  # Gives remainder after "modular arithmetic" ("fmod" in C++)

    # xL = xLag[0] + buff*Ltube 
    # xR = xLag[0] + buff*Ltube + L_AR
    # x = np.logical_and(xPt>=xL,xPt<=xR)
    # af_Val = np.zeros(x.shape)
    
    # af_Val[x] = 0.75*(sin(1/(L_AR)*(2*pi*freq*t+k*xPt[x])-pi/2)+1)


    xL = xLag[0] + buff*Ltube + v*t
    xR = xLag[0] + buff*Ltube + SqWidth + v*t
    xM = (xL+xR)/2
    c = SqWidth/1.5
    x = xPt
    x = np.logical_and(xPt>=xL,xPt<=xR)
    af_Val = np.zeros(x.shape)      
    af_Val[x] = 1                            # Traveling Square Wave
    #af_Val[x] = exp(-(x-xM)**2/(2*c)**2)    # Traveling Gaussian Wave
    
    return af_Val



