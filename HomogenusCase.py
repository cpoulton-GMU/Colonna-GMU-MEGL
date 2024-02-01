import numpy as np
from geometry_tools import hyperbolic
from geometry_tools import drawtools


Degree = 4  # degree of the tree
SegmentLength = np.cos(np.pi/Degree) # branch length
Depth = 6  # number of iterations
Generators = [] # stores the generators for the group associated with the tree. Generators are represented by mobius transforms as 2x2 complex valued matricies
Geodesics = [] # keeps track of edges/geodesics of tree
Points = [] # array of all points/verticies in tree
Origin = np.array([0,1])
ShowPoints = False  # allows you to plot points for verticies
ShowGeodesics = True # allows you to toggle edges on/off
Model = hyperbolic.Model.POINCARE  # uses poincare disk
#Model = hyperbolic.Model.HALFPLANE # for upper half plane
EvenDegree = ((Degree % 2) == 0)  # are we in even or odd case


def FindGenerators():    # this function produces the list of generators as mobius transforms (see page 317 in Flavia's originial paper)
    for i in range(0,Degree):
        if EvenDegree:   # even degree case
            if i == (Degree / 2):
                rot = -1+0j
            else:
                rot = np.clongdouble(np.exp((1j*np.pi*(i))/(Degree / 2))) #lower case zeta in the original paper page 317
            gen = np.array([[1,SegmentLength * rot],[SegmentLength * np.conjugate(rot),1]])
            Generators.append(gen)
        else:    # odd degree case
            rot = np.clongdouble(np.exp((1j * np.pi * (2*i-1)) / (Degree)))
            gen = np.array([[SegmentLength*rot,-1],[1,(-1)*SegmentLength*np.conjugate(rot)]])
            Generators.append(gen)

def FindVerticies(parentwords, parentc, backwardgens, currentDepth):  # function produces array of all verticies from the generators
    newparentwords = [] #new words to be added to the associated group IE compositions of generators
    newparentc = [] #new verticies as points in C
    newbackwardgens = [] #the generator that would take a given parent backwards. Indices need to be identical to newparentwords TODO implement this as dictionary
    for k, parent in enumerate(parentwords):
        parentPoint = hyperbolic.Point([parentc[k].real, parentc[k].imag], model=Model)
        if EvenDegree:    # even degree case
            for i, gen in enumerate(Generators):
                if (i != backwardgens[k]) or (backwardgens[k] == -1): #ignore the generator that would take us backwards toward the origin
                    word = parent@gen #get the next word by composing the current generator with the parent word
                    z = word@Origin #multiply origin by the new word to find the coordinate of the vertex in projective space
                    z = z[0] / z[1] #convert from projective space to C
                    newparentwords.append(word)
                    newbackwardgens.append((i + (Degree/2)) % Degree) #the index of the inverse generator for the generator just applied
                    newparentc.append(z) #save parent as complex number
                    childPoint= hyperbolic.Point([z.real, z.imag], model=Model) #need to convert to special class for geometry_tools
                    Points.append(childPoint)
                    geodesic = hyperbolic.Segment(parentPoint,childPoint)
                    Geodesics.append(geodesic)
        else:
            #TODO ignore generators that take us backwards
            for i, gen in enumerate(Generators):  #odd degree case
                word = parent @ gen
                z = word @ Origin
                z = z[0] / z[1]
                newparentwords.append(word)
                newbackwardgens.append(i)
                newparentc.append(z)
                childPoint = hyperbolic.Point([z.real, z.imag], model=Model)
                Points.append(childPoint)
                geodesic = hyperbolic.Segment(parentPoint, childPoint)
                Geodesics.append(geodesic)


    #base case
    if currentDepth < Depth:
        FindVerticies(newparentwords,newparentc,newbackwardgens,currentDepth+1)


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