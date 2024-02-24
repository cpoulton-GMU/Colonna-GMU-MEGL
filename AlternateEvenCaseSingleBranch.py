import numpy as np
from geometry_tools import hyperbolic
from geometry_tools import drawtools

DEPTH = 100

ODD_DEGREE = 6
EVEN_DEGREE = 4
BOUNDARY_DEGREE = ODD_DEGREE * (EVEN_DEGREE - 1)

ODD_LENGTH = np.cos(np.pi/ODD_DEGREE)
EVEN_LENGTH = np.cos(np.pi/EVEN_DEGREE)
BOUNDARY_LENGTH = 1

MODEL = hyperbolic.Model.POINCARE  # uses poincare disk

Geodesics = []


def Render():    # calls geometry tools to render the hyperbolic tree
    figure = drawtools.HyperbolicDrawing(model=MODEL)
    figure.draw_plane()

    for i, geodesic in enumerate(Geodesics):    # draws the geodesics
        if i == 0:
            geodesic_color = 'RED'
        elif i % 2 == 0:
            geodesic_color = 'GREEN'
        else:
            geodesic_color = 'BLUE'
        figure.draw_geodesic(geodesic, color=geodesic_color)

    figure.show()


def Normalize_Vector(vector):

    return vector[0] / vector[1]

def Add_Geodesic(complex_point_1, complex_point_2):

    complex_point_1 = Normalize_Vector(complex_point_1)
    complex_point_2 = Normalize_Vector(complex_point_2)

    point_1 = hyperbolic.Point([complex_point_1.real, complex_point_1.imag],
                               model=MODEL)  # need to convert to special class for geometry_tools
    point_2 = hyperbolic.Point([complex_point_2.real, complex_point_2.imag],
                               model=MODEL)

    geodesic = hyperbolic.Segment(point_1, point_2)
    Geodesics.append(geodesic)


def Find_Even_Generator():



    upward_direction = np.floor(EVEN_DEGREE / 4)

    r_zeta = np.clongdouble(np.exp((1j * np.pi * upward_direction) / (EVEN_DEGREE / 2)))  # lower case zeta in the original paper page 317
    r_generator = np.array([[1, ODD_LENGTH * r_zeta], [ODD_LENGTH * np.conjugate(r_zeta),
                                                       1]])  # mobius transform for even degree with length of odd degree
    r_point = r_generator @ np.array([0 + 0j, 1 + 0j])
    angle = np.arctan(r_point.real / r_point.imag)

    zeta = np.clongdouble(np.exp(1j * angle[0]))
    generator = np.array([[1, EVEN_LENGTH * zeta], [EVEN_LENGTH * np.conjugate(zeta),
                                                       1]])

    return generator


def Find_Odd_Generator():

    upward_direction = np.floor(ODD_DEGREE / 4)

    zeta = np.clongdouble(np.exp((1j * np.pi * upward_direction) / (ODD_DEGREE / 2)))
    generator = np.array([[1, ODD_LENGTH * zeta], [ODD_LENGTH * np.conjugate(zeta), 1]])

    return generator

def Find_Base_GENERATOR():

    zeta = 0
    generator = np.array([[1, ODD_LENGTH * zeta], [ODD_LENGTH * np.conjugate(zeta), 1]])

    return generator


def Find_Boundry_Geodesic():

    zeta = np.clongdouble(np.exp((1j * np.pi) / (BOUNDARY_DEGREE / 2)))
    generator = np.array([[1, BOUNDARY_LENGTH * zeta], [BOUNDARY_LENGTH * np.conjugate(zeta), 1]])

    boundary_end_point = generator @ np.array([0,1])

    Add_Geodesic(boundary_end_point, np.array([0,1]))

    return


def Calculate_Tree(last_word, even_generator, odd_generator, current_depth, last_point):

    if current_depth % 2 == 0:
        new_word = last_word @ even_generator
    else:
        new_word = last_word @ odd_generator

    new_point = new_word @ np.array([0, 1])
    Add_Geodesic(new_point, last_point)

    if (current_depth < DEPTH):
        current_depth += 1
        Calculate_Tree(new_word, even_generator, odd_generator, current_depth, new_point)


#main
even_generator = Find_Even_Generator()
odd_generator = Find_Odd_Generator()
base_generator = Find_Base_GENERATOR()

Find_Boundry_Geodesic()
Add_Geodesic(np.array([0,1]), base_generator @ np.array([0,1]))

Calculate_Tree(base_generator, even_generator, odd_generator,  0, np.array([ODD_LENGTH, 1 + 0j]))

Render()