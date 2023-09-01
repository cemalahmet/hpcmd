#include "headers/lj_direct_summation.h"

double lj_direct_summation(Atoms &atoms, double epsilon, double sigma) {
    double Ep_tot = 0.0;

    atoms.forces.setZero();

    int n = atoms.nb_atoms();

    for (int i = 0; i < n-1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            // Calculate r, the distance between two atoms
            Positions_t d = atoms.positions.col(j) - atoms.positions.col(i);
            double r = sqrt(d.square().sum());

            // Calculate the lj potential energy and lj force
            double sigma_over_r = sigma / r;

            double sigma_over_r_6 = std::pow(sigma_over_r, 6);
            double sigma_over_r_12 = sigma_over_r_6 * sigma_over_r_6;

            // Total energy of the structure
            double lj_energy = 4 * epsilon * (sigma_over_r_12 - sigma_over_r_6);
            Ep_tot += lj_energy;

            // Update forces on each atom
            double lj_force = -24 * epsilon * (2 * sigma_over_r_12 - sigma_over_r_6) / r;

            atoms.forces.col(i) += d * (lj_force / r);
            atoms.forces.col(j) -= d * (lj_force / r);
        }
    }

    return Ep_tot;
}

double kinetic_energy(const Atoms &atoms, double mass) {
    return mass * atoms.velocities.square().sum() / 2;
}