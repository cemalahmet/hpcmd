#ifndef __LJ_POTENTIAL_H
#define __LJ_POTENTIAL_H

#include "atoms.h"

double lj_direct_summation(Atoms &atoms, double epsilon = 1.0, double sigma = 1.0);
double kinetic_energy(const Atoms &atoms, double mass);

#endif