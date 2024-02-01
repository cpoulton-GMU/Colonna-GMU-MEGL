import numpy as np
from geometry_tools import hyperbolic
from geometry_tools import drawtools


Degree = 4  # degree of the tree
SegmentLength = np.cos(np.pi/Degree) # length of first vertex from origin
Depth = 6  # number of interations
Generators = [] # for even degree, keeps track of free group generators
Geodesics = [] # keeps track of edges/geodesics of tree
Points = [] # array of all points/verticies in tree
Origin = np.array([0,1])
ShowPoints = False  # allows you to plot points for verticies
ShowGeodesics = True # allows you to toggle edges on/off
Model = hyperbolic.Model.POINCARE  # uses poincare disk
#Model = hyperbolic.Model.HALFPLANE # for upper half plane
EvenDegree = ((Degree % 2) == 0)  # boolean that records you are inputting a even degree tree


def FindGenerators():    # this function produces the list of generators
    for i in range(0,Degree):
        if EvenDegree:   # even degree case
            if i == (Degree / 2):
                rot = -1+0j
            else:
                rot = np.clongdouble(np.exp((1j*np.pi*(i))/(Degree / 2)))
            gen = np.array([[1,SegmentLength * rot],[SegmentLength * np.conjugate(rot),1]])
            Generators.append(gen)
        else:    # odd degree case
            rot = np.clongdouble(np.exp((1j * np.pi * (2*i-1)) / (Degree)))
            gen = np.array([[SegmentLength*rot,-1],[1,(-1)*SegmentLength*np.conjugate(rot)]])
            Generators.append(gen)

def FindVerticies(parentwords, parentc, rightmostGenerators, currentDepth):  # function produces array of all verticies from the generators
    newparentwords = []
    newparentc = []
    newrightmostgens = []
    for k, parent in enumerate(parentwords):
        parentPoint = hyperbolic.Point([parentc[k].real, parentc[k].imag], model=Model)
        if EvenDegree:    # even degree case
            for i, gen in enumerate(Generators):
                if (i != ((rightmostGenerators[k] + (Degree / 2)) % Degree)) or (rightmostGenerators[k] == -1):
                    word = parent@gen
                    z = word@Origin
                    z = z[0] / z[1]
                    newparentwords.append(word)
                    newrightmostgens.append(i)
                    newparentc.append(z)
                    childPoint= hyperbolic.Point([z.real, z.imag], model=Model)
                    Points.append(childPoint)
                    geodesic = hyperbolic.Segment(parentPoint,childPoint)
                    Geodesics.append(geodesic)
        else:
            for i, gen in enumerate(Generators):  #odd degree case
                word = parent @ gen
                z = word @ Origin
                z = z[0] / z[1]
                newparentwords.append(word)
                newrightmostgens.append(i)
                newparentc.append(z)
                childPoint = hyperbolic.Point([z.real, z.imag], model=Model)
                Points.append(childPoint)
                geodesic = hyperbolic.Segment(parentPoint, childPoint)
                Geodesics.append(geodesic)



    if currentDepth < Depth:
        FindVerticies(newparentwords,newparentc,newrightmostgens,currentDepth+1)


def Render():    # calls geometry tools to render the hyperbolic tree
    figure = drawtools.HyperbolicDrawing(model=Model)
    figure.draw_plane()

    for geodesic in Geodesics:    # draws the geodesics
        figure.draw_geodesic(geodesic)

    figure.show()


#Main
FindGenerators()
FindVerticies([np.identity(2)],[0],[-1],0)
Render()