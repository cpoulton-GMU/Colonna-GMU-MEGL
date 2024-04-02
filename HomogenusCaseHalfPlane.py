import numpy as np
from geometry_tools import hyperbolic
from geometry_tools import drawtools


DEGREE = 4
TOTAL_ITERATIONS = 3
MODEL = hyperbolic.Model.HALFPLANE
SEGMENT_LENGTH = np.cos(np.pi/DEGREE)
TO_HALFPLANE = np.array([[-1j, 1], [1, -1j]])
ORIGIN = np.array([0, 1])

Generators = []
Geodesics = []


def find_generators():
    for i in range(0, DEGREE):
        angle_of_rot = np.exp((1j*np.pi*i) / (DEGREE / 2))
        generator = np.array([[1, SEGMENT_LENGTH*angle_of_rot], [SEGMENT_LENGTH*np.conjugate(angle_of_rot), 1]])
        Generators.append(generator)


def normalize_projective_vector(proj_coordinate):
    return proj_coordinate / proj_coordinate


def generate_tree(parent_words, iterations):
    ignore = False
    current_words = []
    for parent in parent_words:
        for i in range(0, DEGREE):
            new_word = parent @ Generators[i]
            for potential_match in parent_words:
                print(new_word)
                print(potential_match)
                if new_word == potential_match.all():
                    ignore = True
                    print("ERROR")
                    break
        if not ignore:
            print("TEST")
            current_words.append(new_word)
            child_projective_coordinates = new_word @ ORIGIN
            child_coordinate = normalize_projective_vector(child_projective_coordinates)
            child_point = hyperbolic.Point([child_coordinate.real, child_coordinate.imag], model=MODEL)
            parent_projective_coordinate = parent @ ORIGIN
            parent_coordinate = normalize_projective_vector(parent_projective_coordinate)
            parent_point = hyperbolic.Point([parent_coordinate.real, parent_coordinate.imag], model=MODEL)
            geodesic = hyperbolic.Segment(parent_point, child_point)
            Geodesics.append(geodesic)
        else:
            ignore = False

    if iterations < TOTAL_ITERATIONS:
        generate_tree(current_words, iterations + 1)


def render():    # calls geometry tools to render the hyperbolic tree
    figure = drawtools.HyperbolicDrawing(model=MODEL)
    figure.draw_plane()

    for geodesic in Geodesics:    # draws the geodesics
        print(geodesic)
        figure.draw_geodesic(geodesic)

    figure.show()


#MAIN
find_generators()
generate_tree([TO_HALFPLANE], 0)
render()








