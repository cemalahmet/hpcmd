#include "headers/atoms.h"
#include "headers/berendsen_thermostat.h"
#include "headers/ducastelle.h"
#include "headers/lj_cutoff.h"
#include "headers/lj_direct_summation.h"
#include "headers/verlet.h"
#include "headers/xyz.h"

#include <cmath>
#include <iostream>
#include <sys/time.h>

// IMPORTANT: need to call ./mackay/create_objects.sh script
// before running this
int main(int argc, char *argv[]) {
    // Constants & settings
    constexpr double mass = 20405.736652;
    constexpr double cutoff = 5.0;
    constexpr double kB = 8.617333262e-5;
    constexpr double delta_q = 25;
    constexpr int t_relax = 1000;

    double total_time = 40000;
    double timestep = 1;
    long nb_steps = floor(total_time / timestep);

    // Get data for energy vs. temperature on cluster_923
    std::ofstream traj("traj.xyz");
    std::ofstream out("out.dat");

    auto [names, positions]{read_xyz("cluster_923.xyz")};
    Atoms atoms{positions};

    NeighborList neighbor_list;

    double temp_avg = 0.0;
    double total_energy_avg = 0.0;
    for (int i = 0; i < nb_steps; i++) {
        verlet_step1(atoms, mass, timestep);

        neighbor_list.update(atoms, cutoff);
        double ep = ducastelle(atoms, neighbor_list, cutoff);

        verlet_step2(atoms, mass, timestep);

        // Measure for t_relax when another t_relax passed after heat injection
        if (i % (2 * t_relax) >= t_relax) {
            double ek = kinetic_energy(atoms, mass);
            double temp = temperature(ek, atoms.nb_atoms(), kB);

            temp_avg += temp;
            total_energy_avg += ek + ep;
        }

        // Compute the data and inject heat every 2t_relax time)
        if (i && i % (2 * t_relax) == 0) {
            double ek = kinetic_energy(atoms, mass);

            // compute average temp
            temp_avg /= t_relax;
            total_energy_avg /= t_relax;

            // for the plot of energy vs. temperature
            out << temp_avg << std::endl;
            out << total_energy_avg << std::endl;

            // add energy of delta Q (rescale velocities)
            double lambda = sqrt((delta_q / ek) + 1);
            atoms.velocities = atoms.velocities * lambda;

            temp_avg = 0.0;
            total_energy_avg = 0.0;
        }

        //  show trajectory
        if (i % 10 == 0) {
            write_xyz(traj, atoms);
            if (i % 1000 == 0)
                std::cout << i / 1000 << std::endl;
        }
    }

    out.close();
    traj.close();

    // Now output the heat capacity latent heat for
    // different cluster sizes. I have not found a way to automatically find the
    // melting point. That's why I use the data from here and analyze it myself.
    // Make sure to run this before:
    // $ sh ../../mackay/create_objects.sh
    double delta_qs[13] = {0, 0, 1, 1, 5, 10, 20, 40, 60, 80, 110, 120, 150};

    for (int i = 3; i <= 12; i++) {
        double delta_q2 = delta_qs[i];

        std::ofstream out2("out" + std::to_string(i) + ".dat");
        std::string filename = "../../mackay/" + std::to_string(i) + ".xyz";

        auto [names2, positions2] = read_xyz(filename);
        Atoms atoms2{positions2};

        int N = atoms2.nb_atoms();
        out2 << N << std::endl;

        NeighborList neighbor_list2;

        double total_energy_avg2 = 0.0;
        double temp_avg2 = 0.0;

        double prev_temp_avg = 0.0;
        double prev_energy_avg = 0.0;

        int t_relax2 = 1000;

        for (int j = 0; j < 50000; j++) {
            verlet_step1(atoms2, mass, timestep);

            neighbor_list2.update(atoms2, cutoff);
            double ep = ducastelle(atoms2, neighbor_list2, cutoff);

            verlet_step2(atoms2, mass, timestep);

            // Measure for t_relax when another t_relax passed after heat
            // injection
            if (j % (2 * t_relax2) >= t_relax2) {
                double ek = kinetic_energy(atoms2, mass);
                double temp = temperature(ek, atoms2.nb_atoms(), kB);

                temp_avg2 += temp;
                total_energy_avg2 += ek + ep;
            }

            // Compute the data and inject heat every 2t_relax time)
            if (j && j % (2 * t_relax2) == 0) {
                double ek = kinetic_energy(atoms2, mass);
                // compute average temp and energy
                temp_avg2 /= t_relax2;
                total_energy_avg2 /= t_relax2;

                if (prev_temp_avg != 0) {
                    double heat_cap = (total_energy_avg2 - prev_energy_avg) /
                                      (temp_avg2 - prev_temp_avg);

                    out2 << ((j + 0.0) / (2 * t_relax2)) * delta_q2 << std::endl; // heat injected so far
                    out2 << temp_avg2 << std::endl;  // temperature
                    out2 << heat_cap << std::endl; // heat capacity
                    out2 << std::endl;
                }

                prev_energy_avg = total_energy_avg2;
                prev_temp_avg = temp_avg2;

                temp_avg2 = 0.0;
                total_energy_avg2 = 0.0;

                // add energy of delta Q (rescale velocities)
                double lambda = sqrt((delta_q2 / ek) + 1);
                atoms2.velocities = atoms2.velocities * lambda;

            }
        }
        out2.close();
    }
}
