import ctypes
import os

import numpy as np

lib_path = os.path.abspath("libtsp.so")
tsp = ctypes.CDLL(lib_path)

tsp.farthestInsertion.argtypes = (
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

tsp.nearestInsertion.argtypes = (
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

tsp.localSearch2opt.argtypes = (
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int),  # cost matrix
    ctypes.POINTER(ctypes.c_int),  # tour (initial tour used also for output)
)
def nearestInsertion(cost_matrix):
    cost_matrix = np.array(cost_matrix, dtype=np.int32)
    size = cost_matrix.shape[0]
    cost_mat_flat = np.ascontiguousarray(cost_matrix)
    tour = np.zeros(size + 1, dtype=np.int32)

    mat_ptr = cost_mat_flat.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    tour_ptr = tour.ctypes.data_as(ctypes.POINTER(ctypes.c_int))

    tsp.nearestInsertion(ctypes.c_int(size), mat_ptr, tour_ptr)
    return tour.tolist()

def farthestInsertion(cost_matrix):
    cost_matrix = np.array(cost_matrix, dtype=np.int32)
    size = cost_matrix.shape[0]
    cost_mat_flat = np.ascontiguousarray(cost_matrix)
    tour = np.zeros(size + 1, dtype=np.int32)

    mat_ptr = cost_mat_flat.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    tour_ptr = tour.ctypes.data_as(ctypes.POINTER(ctypes.c_int))

    tsp.farthestInsertion(ctypes.c_int(size), mat_ptr, tour_ptr)
    return tour.tolist()


def localSearch(cost_matrix, tour):
    cost_matrix = np.array(cost_matrix, dtype=np.int32)
    size = cost_matrix.shape[0]
    cost_mat_flat = np.ascontiguousarray(cost_matrix)
    tour = np.array(tour, dtype=np.int32)

    mat_ptr = cost_mat_flat.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    tour_ptr = tour.ctypes.data_as(ctypes.POINTER(ctypes.c_int))

    tsp.localSearch2opt(ctypes.c_int(size), mat_ptr, tour_ptr)
    return tour.tolist()


def solve_heuristic(cost_matrix):
    tour = nearestInsertion(cost_matrix)
    tour = localSearch(cost_matrix, tour)
    return tour
