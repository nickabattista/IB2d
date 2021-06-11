from math import cos, sin, pi, sqrt, exp
import numpy as np
import matplotlib.pyplot as plt

def test_Activation_Wave():

    xLag = np.arange(1,3.00001,0.025)


    freq = 10                      # frequency of traveling wave down tube
    Ltube = xLag[-1] - xLag[0]     # length of heart tube
    buff = 0.25                    # percent-buff on each end of heart tube
    L_AR = (1-2*buff)*Ltube        # length of the actual activation region of tube
    SqWidth = L_AR/10              # width of the square traveling activation wave 
    v = (L_AR-SqWidth) * freq      # traveling wave velocity
    k = 2*pi*freq/v                # Wave-number for traveling wave

    
    time = np.arange(0,1.0001,0.001)

    for jj in range(time.size):
        
        t = time[jj]
        t = np.fmod(t,1/freq) # gives remainder after "modular arithmetic" ("fmod" in C++)

        xL = xLag[0]+ buff*Ltube + v*t
        xR = xLag[0]+ buff*Ltube + SqWidth + v*t
        xM = (xL+xR)/2
        c = SqWidth/1.5
        
        af_Val = np.zeros(len(xLag))
        for ii in range(len(xLag)):
            x = xLag[ii]
            if x >= xL and x <= xR:
                af_Val[ii] = 1                        # Traveling Square Wave
                #af_Val = exp(-(x-xM)**2/(2*c**2))    # Traveling Gaussian Wave
            else:
                af_Val[ii] = 0.0
            
        
        plt.plot(xLag,af_Val,'r-')
        plt.plot(xLag,af_Val,'*')
        strTitle = 'Time(s) = {}'.format(round(t,3))
        plt.title(strTitle)
        plt.axis([0.8, 3.2, -0.2, 1.2])
        plt.draw()
        plt.pause(0.1)
        plt.clf

if __name__ == "__main__":
    test_Activation_Wave()