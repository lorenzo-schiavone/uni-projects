#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <iostream>

void apply2optMove(int* tour, const std::pair<int, int>& move) {
    int a = move.first;
    int b = move.second;
    while (a < b) {
            std::swap(tour[a], tour[b]);
            a++;
            b--;
    }
}

inline int &A_at(int *A_flat, int nrows, int i, int j){
    return A_flat[i * nrows + j];
}

int findBestNeighbor(int size, int* cost_matrix, int* tour, std::pair<int, int>& move) {
    int bestCostVariation = 1000000;
    for (int a = 1; a < size - 1; a++) {
        int h = tour[a - 1];
        int i = tour[a];
        for (int b = a + 1; b < size; b++) {
            int j = tour[b];
            int l = tour[b + 1];
            int neighCostVariation = - A_at(cost_matrix, size, h, i) - A_at(cost_matrix, size, j, l) + A_at(cost_matrix, size, h, j) + A_at(cost_matrix, size, i, l);
            if (neighCostVariation < bestCostVariation) {
                bestCostVariation = neighCostVariation;
                move.first = a;
                move.second = b;
            }
        }
    }
    return bestCostVariation;
}
extern "C" {
void localSearch2opt(int size, int* cost_matrix, int* tour) {
    std::pair<int, int> move;
    while (true) {
        int improvement2opt = findBestNeighbor(size, cost_matrix, tour, move);
        if (improvement2opt < 0) {
            apply2optMove(tour, move);
        } else {
            break;
        }
    }
}


void farthestInsertion(int size, int* cost_matrix, int* tour){

    std::vector<bool> unvisited(size, true);
    unvisited[0]=false;
    int i_best = -1;
    int max_dist = -100000;

    for ( int i = 1 ; i < size-1; ++i ) {
            if (A_at(cost_matrix, size, 0, i) > max_dist) {
                max_dist = A_at(cost_matrix, size, 0, i);
                i_best = i;
            }
    }
    unvisited[i_best] = false;
    std::vector<int> local_tour;
    local_tour.push_back(0);
    local_tour.push_back(i_best);
    local_tour.push_back(0);

    int local_size = 3;

    while (local_size<size+1) {
        int r = -1;
        int best_value = -1000000;

        for (int v=1; v<size; v++) {
            if (!unvisited[v])continue;
            int max_dist_to_tour = 0;
            for (int j=0; j<local_size-1; j++) {
                int w=local_tour[j];
                max_dist_to_tour = std::max(max_dist_to_tour, A_at(cost_matrix, size, v, w));
            }
            if (max_dist_to_tour > best_value) {
                best_value = max_dist_to_tour;
                r = v;
            }
        }

        int best_pos = -1;
        int best_cost = 1000000;

        for (int idx = 0; idx < local_size - 1; idx++) {
            int i = local_tour[idx], j = local_tour[idx + 1];
            int insertion_cost = A_at(cost_matrix, size, i, r) + A_at(cost_matrix, size, r, j) - A_at(cost_matrix, size, i, j);

            if (insertion_cost < best_cost) {
                best_cost = insertion_cost;
                best_pos = idx + 1;
            }
        }
        local_tour.insert(local_tour.begin() + best_pos, r);
        local_size++;
        unvisited[r]=false;
    }
    for (int i = 0; i < size+1; ++i) {
            tour[i] = local_tour[i];
    }
}

void nearestInsertion(int size, int* cost_matrix, int* tour) {
    std::vector<bool> unvisited(size, true);
    unvisited[0] = false;

    int i_best = -1;
    int min_dist = 1000000;

    for (int i = 1; i < size; ++i) {
        int dist = A_at(cost_matrix, size, 0, i);
        if (dist < min_dist) {
            min_dist = dist;
            i_best = i;
        }
    }

    unvisited[i_best] = false;
    std::vector<int> local_tour;
    local_tour.push_back(0);
    local_tour.push_back(i_best);
    local_tour.push_back(0);
    int local_size = 3;

    while (local_size < size + 1) {
        int r = -1;
        int best_selection_value = 1000000;

        for (int v = 1; v < size; v++) {
            if (!unvisited[v]) continue;

            int min_dist_to_tour = 1000000;
            for (int idx = 0; idx < local_size - 1; idx++) {
                int w = local_tour[idx];
                int dist = A_at(cost_matrix, size, v, w);
                if (dist < min_dist_to_tour) {
                    min_dist_to_tour = dist;
                }
            }

            if (min_dist_to_tour < best_selection_value) {
                best_selection_value = min_dist_to_tour;
                r = v;
            }
        }

        int best_pos = -1;
        int best_insertion_cost = 1000000;

        for (int idx = 0; idx < local_size - 1; idx++) {
            int i = local_tour[idx];
            int j = local_tour[idx + 1];

            int insertion_cost = A_at(cost_matrix, size, i, r) + A_at(cost_matrix, size, r, j) - A_at(cost_matrix, size, i, j);

            if (insertion_cost < best_insertion_cost) {
                best_insertion_cost = insertion_cost;
                best_pos = idx + 1;
            }
        }

        local_tour.insert(local_tour.begin() + best_pos, r);
        local_size++;
        unvisited[r] = false;
    }

    for (int i = 0; i < size + 1; ++i) {
        tour[i] = local_tour[i];
    }
}
}
