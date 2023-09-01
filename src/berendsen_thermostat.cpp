#include "headers/berendsen_thermostat.h"
#include "headers/lj_direct_summation.h"

#include <cmath>
#include <iostream>

// Return the temperature given the atoms object
double temperature(Atoms &atoms, double mass, double kB) {
    return temperature(kinetic_energy(atoms, mass), atoms.nb_atoms(), kB);
}

// Return the temperature given the kinetic energy and number of atoms
double temperature(double Ek, size_t n, double kB) {
    return Ek / (double(n) * 1.5 * kB);
}

void berendsen_thermostat(Atoms &atoms, double mass, double goal_temperature, double timestep,
                          double relaxation_time, double kB) {
    double T = temperature(atoms, mass, kB);

    // To avoid dividing by 0
    if (T == 0) return;

    // Rescale all velocities
    double lambda = sqrt(1 + ((goal_temperature / T) - 1) * (timestep / relaxation_time));
    atoms.velocities *= lambda;
}


