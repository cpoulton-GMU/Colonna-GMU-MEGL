import numpy as np
from geometry_tools import hyperbolic
from geometry_tools import drawtools


ODD_DEGREE = 6  # odd alternating degree.
EVEN_DEGREE = 4 # even alternating degree
BOUNDRY_DEGREE = 18 # boundry tree
SEGMENTLENGTH = np.cos(np.pi/BOUNDRY_DEGREE) # branch length
DEPTH = 1 # number of iterations
Generators_even = dict()- # stores the generators for the group associated with the tree. Generators are represented by mobius transforms as 2x2 complex valued matricies
Generators_odd = dict()
Generators_boundry = dict()
AlternatingGeodesics = [] # keeps track of edges/geodesics of tree
BoundryGeodesics = []
Points = [] # array of all points/verticies in tree
Origin = np.array([0,1])
ShowPoints = False  # allows you to plot points for verticies
ShowGeodesics = True # allows you to toggle edges on/off
Model = hyperbolic.Model.POINCARE  # uses poincare disk
#Model = hyperbolic.Model.HALFPLANE # for upper half plane
EvenDegree = True


def FindGenerators(Generators, Degree):    # this function produces the list of generators as mobius transforms (see page 317 in Flavia's originial paper)
    for i in range(0,Degree):
        if EvenDegree:   # even degree case
            if i == (Degree / 2):
                rot = -1+0j
            else:
                rot = np.clongdouble(np.exp((1j*np.pi*(i))/(Degree / 2))) #lower case zeta in the original paper page 317
            gen = np.array([[1,SEGMENTLENGTH * rot],[SEGMENTLENGTH * np.conjugate(rot),1]])
            Generators.update({i : gen})
        else:    # odd degree case
            rot = np.clongdouble(np.exp((1j * np.pi * (2*i-1)) / (Degree)))
            gen = np.array([[SEGMENTLENGTH*rot,-1],[1,(-1)*SEGMENTLENGTH*np.conjugate(rot)]])
            Generators.update({i : gen})
        
    return Generators

def Norm_Vector(vector):
    return vector[0] / vector[1]

def Find_Gen_Inverse_Even(generator_index, total_gens):
    return (generator_index + (total_gens/2)) % total_gens

def FindVerticies(Generators, Parents, Calculate_Inverse, CurrentDepth):  # function produces array of all verticies from the generators
    newparents = [] #[[word, complex coordinate, backward generator]]
    for parent in Parents:
        print(parent[1])
        parentPoint = hyperbolic.Point([parent[1].real, parent[1].imag], model=Model)
        if EvenDegree:    # even degree case
            for i in range(len(Generators)):
                gen = Generators.get(i)
                new_inverse = Find_Gen_Inverse_Even(i,len(Generators))
                if (i != parent[2]) or (Calculate_Inverse == False): #ignore the generator that would take us backwards toward the origin
                    word = parent[0]@gen #get the next word by composing the current generator with the parent word
                    z = word@Origin #multiply origin by the new word to find the coordinate of the vertex in projective space
                    z = Norm_Vector(z)
                    newparents.append([word,z,new_inverse])
                    childPoint= hyperbolic.Point([z.real, z.imag], model=Model) #need to convert to special class for geometry_tools
                    Points.append(childPoint)
                    geodesic = hyperbolic.Segment(parentPoint,childPoint)
                    BoundryGeodesics.append(geodesic)

    #base case
    if CurrentDepth < DEPTH:
        print("test2")
        FindVerticies(Generators, newparents, True, CurrentDepth + 1)


def FindVerticiesAlternating(Generators, Parents, Calculate_Inverse, CurrentDepth):  # function produces array of all verticies from the generators
    newparents = [] #[[word, complex coordinate, backward generator]]
    for parent in Parents:
        print(parent[1])
        parentPoint = hyperbolic.Point([parent[1].real, parent[1].imag], model=Model)
        if EvenDegree:    # even degree case
            for i in range(len(Generators)):
                gen = Generators.get(i)
                new_inverse = Find_Gen_Inverse_Even(i,len(Generators))
                if (i != parent[2]) or (Calculate_Inverse == False): #ignore the generator that would take us backwards toward the origin
                    word = parent[0]@gen #get the next word by composing the current generator with the parent word
                    z = word@Origin #multiply origin by the new word to find the coordinate of the vertex in projective space
                    z = Norm_Vector(z)
                    newparents.append([word,z,new_inverse])
                    childPoint= hyperbolic.Point([z.real, z.imag], model=Model) #need to convert to special class for geometry_tools
                    Points.append(childPoint)
                    geodesic = hyperbolic.Segment(parentPoint,childPoint)
                    AlternatingGeodesics.append(geodesic)


    #base case
    if CurrentDepth < DEPTH:
        if (CurrentDepth % 2) == 0:
            newgens = Generators_even
        else:
            newgens = Generators_odd
        FindVerticiesAlternating(newgens, newparents, True, CurrentDepth + 1)


def Render(Color, Geodesics):    # calls geometry tools to render the hyperbolic tree
    figure = drawtools.HyperbolicDrawing(model=Model)
    figure.draw_plane()

    for geodesic in Geodesics:    # draws the geodesics
        figure.draw_geodesic(geodesic, color=Color)

    figure.show()


#Main
Generators_even = FindGenerators(dict(), EVEN_DEGREE)
Generators_odd = FindGenerators(dict(), ODD_DEGREE)
Generators_boundry = FindGenerators(dict(), BOUNDRY_DEGREE)


InitialParentAlt =  [[Generators_odd[0], Norm_Vector(Generators_odd[0]@np.array([0+0j,1+0j])), Find_Gen_Inverse_Even(0,len(Generators_odd))]]
FindVerticiesAlternating(Generators_odd, InitialParentAlt, True, 0)
InitialParentBoundry =  [[Generators_boundry[0], Norm_Vector(Generators_boundry[0]@np.array([0+0j,1+0j])), Find_Gen_Inverse_Even(0,len(Generators_boundry))]]
FindVerticies(Generators_boundry, InitialParentBoundry, True, 0)

figure = drawtools.HyperbolicDrawing(model=Model)
figure.draw_plane()

for geodesic in BoundryGeodesics:
    figure.draw_geodesic(geodesic, color='RED')

for geodesic in AlternatingGeodesics:
    figure.draw_geodesic(geodesic, color='BLACK')





figure.show()