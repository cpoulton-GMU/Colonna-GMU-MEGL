import sys
import numpy as np
from geometry_tools import hyperbolic
from geometry_tools import drawtools


DEGREE = int(sys.argv[1])
TOTAL_ITERATIONS = int(sys.argv[2])
MODEL = hyperbolic.Model.POINCARE
SEGMENT_LENGTH = np.cos(np.pi/DEGREE)
ORIGIN = np.array([0, 1])

Generators = []
Geodesics = []


class VertexWord:
    def __init__(self, projection_matrix, backward_projection_matrix):
        self.proj_mtrx = projection_matrix #matrix representation of the disk automorphisim associated with the vertex
        self.bk_proj_mtrx = backward_projection_matrix #matrix associated with the inverse generator of this word's rightmost generator (IE the generator that takes us backwards)


def find_generators():#find the matrix associated with each generator of the free group associated with the given degree, then append it's inverse generator.
    generators = []
    for i in range(0, DEGREE):
        angle_of_rot = np.exp((1j*np.pi*i) / (DEGREE / 2))
        generator = np.array([[1, SEGMENT_LENGTH*angle_of_rot], [SEGMENT_LENGTH*np.conjugate(angle_of_rot), 1]]) #see Colona's paper on hyperbolic tree embeddings, page 317 for the derivation of this formula
        generators.append(generator)

    for i, gens in enumerate(generators):
        inverse_id = ((i + (DEGREE / 2)) % DEGREE)
        inverse_id = int(inverse_id)
        Generators.append(VertexWord(gens, generators[inverse_id]))


def normalize_projective_vector(proj_coordinate): #convert from points in CP1 to points in the complex plane
    return proj_coordinate[0] / proj_coordinate[1]


def generate_tree(parent_words, iterations): #iterative right actions by each generator
    current_words = []
    for parent in parent_words:
        for gens in Generators:
            if not np.array_equal(gens.bk_proj_mtrx, parent.bk_proj_mtrx): #ignore the generator that sends us backwards
                new_word = parent.proj_mtrx @ gens.proj_mtrx
                current_words.append(VertexWord(new_word, gens.bk_proj_mtrx))
                child_projective_coordinates = new_word @ ORIGIN
                child_coordinate = normalize_projective_vector(child_projective_coordinates)
                child_point = hyperbolic.Point([child_coordinate.real, child_coordinate.imag], model=MODEL)
                parent_projective_coordinate = parent.proj_mtrx @ ORIGIN
                parent_coordinate = normalize_projective_vector(parent_projective_coordinate)
                parent_point = hyperbolic.Point([parent_coordinate.real, parent_coordinate.imag], model=MODEL)
                geodesic = hyperbolic.Segment(parent_point, child_point)
                Geodesics.append(geodesic)

    if iterations < TOTAL_ITERATIONS:
        generate_tree(current_words, iterations + 1)


def render():    # calls geometry tools to render the hyperbolic tree
    figure = drawtools.HyperbolicDrawing(model=MODEL)
    figure.draw_plane()

    for geodesic in Geodesics:    # draws the geodesics
        figure.draw_geodesic(geodesic)

    figure.show()


#MAIN
find_generators()
generate_tree([VertexWord(np.identity(2), np.identity(2))], 0)
render()