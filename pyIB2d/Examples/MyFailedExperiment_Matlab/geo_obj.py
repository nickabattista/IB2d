# -*- coding: utf-8 -*-
# Python Class to create geometry objects for
# use with IB2D and IBAMR
# ==============================================================================
# Classes to handle the various geo_objects
#      - Vertex: the x,y grid of points. It's elements will be used as
#           as references for other objects in the plane.
#      - Spring
#      - Beam
#      - Porous Point
#      - Mass
# ==============================================================================

class Vertex():
    """
    defines the x,y Lagrangian nodes
    """
    def __init__(self,x,y):
        self.x = x
        self.y = y

    def __eq__(self,other):
        return (self.x == other.x) and (self.y == other.y)

    def __lt__(self,other):
        if self.x==other.x:
            if self.y<other.y:
                return True
            else:
                return False
        elif self.x<self.y:
            return True
        else:
            return False

    def getPos(self):
        return (self.x, self.y)

    def getType(self):
        return "vertex"

    def printString(self):
        """
        Print vertex string for the .node file.
        """
        return repr(self.x) + " " + repr(self.y) + "\n" # neede to implement

class Spring():
    """
    defines a spring element.
    """
    def __init__(self,mastId,slaveId, stiff=2.5e4, restlen=0):
        self.master = mastId    # master node ID
        self.slave = slaveId    # slave node ID
        self.stiff = stiff      # spring stiffness
        self.restlen = restlen  # resting length
        #self.beta = beta        # deg. of nonlinearity

    def getType(self):
        return "spring"

    def printString(self):
        spring_str = repr(self.master) + " " + repr(self.slave) + \
                    " " + repr(self.stiff) + " " + repr(self.restlen) + "\n"
        return spring_str

class Beam():
    """
    defines a beam element.
    """
    def __init__(self, lID, mID, rID, stiff, curv):
        self.lID = lID      # left node in beam
        self.mID = mID      # middle node in beam
        self.rID = rID      # right node in beam
        self.kb = stiff     # beam stiffness
        self.c  = curv      # beam curvature

    def getType(self):
        return "beam"

    def printString(self):
        beam_str = repr(self.lID) + " " + repr(self.mID) + " " + \
                   repr(self.rID) + " " + repr(self.kb)  + " " + \
                   repr(self.c) + "\n"
        return beam_str

class PorousPt():
    """
    defines a porous point.
    ONLY works in IB2d, does NOT work in IBAMR.
    """
    def __init__(self, lagID, pcoeff, stenID):
        self.lagID  = lagID     # node which is porous
        self.pcoeff = pcoeff    # porosity coefficient
        self.stenID = stenID    # stencil ID

    def getType(self):
        return "porous"

    def printString(self):
        pstring = repr(self.lagID) + " " + repr(self.pcoeff) + " " + \
                  repr(self.stenID) + "\n"
        return pstring

class Mass():
    """
    defines a mass.
    """
    def __init__(self, lagID, stiff, kg):
        self.lagID = lagID  # Lagrangian node ID
        self.stiff = stiff  # stiffness of Mass
        self.kg    = kg     # mass of point

    def getType(self):
        return "mass"

    def printString(self):
        mstring = repr(self.lagID) + " " + repr(self.stiff) + " " + \
               repr(self.kg) + "\n"
        return mstring

#===============================================================================
# Function to write the various geometry files
#===============================================================================
def writeFile(filename, geo_list):
    """
    writes the .OBJ files based on filename.
    """
    fname = filename + "." + geo_list[0].getType()
    f = open(fname, "w")
    n = len(geo_list)
    f.write(repr(n)+"\n")
    for elem in geo_list:
        f.write(elem.printString())
    f.close()

#===============================================================================
# Testing function to make sure that this file is working correctly.
#===============================================================================
def test():
    print("Beginning test.")
    v1 = Vertex(3,4)
    v2 = Vertex(5,6)
    vlist = [v1, v2]
    writeFile('test', vlist)

    s1 = Spring(3,4, 1e-4, 1.1, 4)
    s2 = Spring(5,6, 1e-4, 1.1, 4)
    slist = [s1, s2]
    writeFile('test', slist)
    print("Test completed.")

if __name__ == '__main__': test()
