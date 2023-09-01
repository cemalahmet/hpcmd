#include "headers/verlet.h"

void verlet_step1(Positions_t &positions, Velocities_t &velocities,
                  const Velocities_t &forces, double mass, double timestep) {
    // Update Velocity
    velocities += forces * (0.5 * timestep / mass);
    // New location
    positions += velocities * timestep;
}

void verlet_step2(Velocities_t &velocities, const Velocities_t &forces, double mass, double timestep) {
    // Update Velocity
    velocities += forces * (0.5 * timestep / mass);
}


void verlet_step1(Atoms &atoms, double mass, double timestep) {
    verlet_step1(atoms.positions, atoms.velocities, atoms.forces, mass, timestep);
}
void verlet_step2(Atoms &atoms, double mass, double timestep) {
    verlet_step2(atoms.velocities, atoms.forces, mass, timestep);
}


