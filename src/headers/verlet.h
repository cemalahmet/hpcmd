#ifndef __VERLET_H
#define __VERLET_H

#include "atoms.h"
#include "types.h"
#include <Eigen/Dense>

void verlet_step1(Positions_t &positions, Velocities_t &velocities,
                  const Velocities_t &forces, double mass, double timestep);
void verlet_step2(Velocities_t &velocities, const Velocities_t &forces,
                  double mass, double timestep);

void verlet_step1(Atoms &atoms, double mass, double timestep);
void verlet_step2(Atoms &atoms, double mass, double timestep);

#endif // __VERLET_H