//
// Created by ahmet on 19/08/23.
//

#ifndef __LJ_CUTOFF_H
#define __LJ_CUTOFF_H


#include "atoms.h"
#include "neighbors.h"


double lj_cutoff(Atoms &atoms, NeighborList &neighbor_list, double epsilon, double sigma, double cutoff);

#endif // __LJ_CUTOFF_H
