#include "headers/atoms.h"
#include "headers/lj_direct_summation.h"
#include "headers/verlet.h"
#include "headers/xyz.h"

#include <iostream>
#include <fstream>
#include <cmath>

int main(int argc, char *argv[]) {
    auto [names, positions, velocities]{read_xyz_with_velocities("lj54.xyz")};
    Atoms atoms{positions, velocities};

    // See the conservation of energy
    std::ofstream traj("traj.xyz");
    std::ofstream out("out.dat");

    if(!out  || out.fail() || !traj || traj.fail()) {
        std::cerr << "Error: file could not be opened" << std::endl;
        return -1;
    }

    constexpr double sigma = 1.0;
    constexpr double epsilon = 1.0;
    constexpr double mass = 1.0;

    double total_time = 100 * sqrt(mass * sigma * sigma / epsilon);
    double timestep = 0.001 * sqrt(mass * sigma * sigma / epsilon);
    long nb_steps = floor(total_time / timestep);

    for (int i = 0; i < nb_steps; i++) {
        verlet_step1(atoms, mass, timestep);

        double Ep = lj_direct_summation(atoms, epsilon, sigma);

        verlet_step2(atoms, mass, timestep);
        double Ek = kinetic_energy(atoms, mass);

        write_xyz(traj, atoms);

        // Output the timestep index, potential energy, and kinetic energy
        out << i << std::endl;
        out << Ep << std::endl;
        out << Ek << std::endl;
    }
    traj.close();
    out.close();

    // Test total energy vs timesteps
    std::ofstream out2("out2.dat");

    int ts_count = 6;
    double ts_list[6] = {0.001, 0.002, 0.005, 0.01, 0.02, 0.03};

    for (int k = 0; k < ts_count; k++) {
        double ts = ts_list[k];
        long nbs = floor(total_time / ts);

        // Re-read the file
        auto [names2, positions2, velocities2]{read_xyz_with_velocities("lj54.xyz")};
        Atoms atoms2{positions2, velocities2};

        for (int i = 0; i < nbs; i++) {
            verlet_step1(atoms2, mass, ts);

            double Ep = lj_direct_summation(atoms2, epsilon, sigma);

            verlet_step2(atoms2, mass, ts);
            double Ek = kinetic_energy(atoms2, mass);

            // output total energy for each timestep
            out2 << Ep + Ek << std::endl;
        }
        // indicator for changing to the next timestep value
        out2 << "###" << std::endl;
    }
    out2.close();
}
