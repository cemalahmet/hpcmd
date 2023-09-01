#include "headers/atoms.h"
#include "headers/lj_direct_summation.h"
#include "headers/verlet.h"
#include "headers/xyz.h"
#include "headers/berendsen_thermostat.h"

#include <cmath>
#include <iostream>
#include <sys/time.h>

int main(int argc, char *argv[]) {
    constexpr double sigma = 1.0;
    constexpr double epsilon = 1.0;
    constexpr double mass = 1.0;

    double lattice_constant{1.4 * sigma};

    double total_time = 100 * sigma * sqrt(mass / epsilon);
    double timestep = 0.002 * sigma * sqrt(mass / epsilon);

    double goal_temperature = 0.35;
    double kB = 1.0;

    long nb_steps = floor(total_time / timestep);

    std::ofstream traj("traj.xyz");

    double relaxation_times[3]{0.02, 0.1, 0.5};

    for (int l = 0; l < 3; l++) {
        double relaxation_time = relaxation_times[l];

        std::ofstream out2("out" + std::to_string(l + 1) + ".dat");
        Atoms atoms2{size_t(125)};

        // Cubic lattice
        for (int x{0}, i{0}; x < 5; ++x) {
            for (int y{0}; y < 5; ++y) {
                for (int z{0}; z < 5; ++z, ++i) {
                    atoms2.positions(0, i) += x * lattice_constant;
                    atoms2.positions(1, i) += y * lattice_constant;
                    atoms2.positions(2, i) += z * lattice_constant;
                }
            }
        }

        for (int i = 0; i < nb_steps; i++) {
            verlet_step1(atoms2, mass, timestep);

            lj_direct_summation(atoms2, epsilon, sigma);

            verlet_step2(atoms2, mass, timestep);

            berendsen_thermostat(atoms2, mass, goal_temperature, timestep,relaxation_time, kB);

            if (l == 0) write_xyz(traj, atoms2);
            out2 << temperature(atoms2, mass, kB) << std::endl;
        }

        traj.close();
        out2.close();

    }

    std::ofstream out("out.dat");

    struct timeval start, end;
    for (int N = 2; N <= 6; N++) {
        for (int M = 0; M < 3; M++) {
            int N1 = N;
            int N2 = (N + (M > 0));
            int N3 = (N + (M > 1));

            Atoms atoms{size_t(N1 * N2 * N3)};

            // Cubic lattice
            for (int x{0}, i{0}; x < N1; ++x) {
                for (int y{0}; y < N2; ++y) {
                    for (int z{0}; z < N3; ++z, ++i) {
                        atoms.positions(0, i) += x * lattice_constant;
                        atoms.positions(1, i) += y * lattice_constant;
                        atoms.positions(2, i) += z * lattice_constant;
                    }
                }
            }

            gettimeofday(&start, NULL);
            for (int i = 0; i < nb_steps; i++) {
                verlet_step1(atoms, mass, timestep);

                lj_direct_summation(atoms, epsilon, sigma);

                verlet_step2(atoms, mass, timestep);

                if (i % 10 == 0)
                    berendsen_thermostat(atoms, mass, goal_temperature, timestep,10*timestep, kB);
            }
            gettimeofday(&end, NULL);
            out << N1 * N2 * N3 << std::endl;
            std::cout << N1 * N2 * N3 << " "
                      << (end.tv_sec - start.tv_sec) * 1000.0 +
                             (end.tv_usec - start.tv_usec) / 1000.0
                      << std::endl;
            out << (end.tv_sec - start.tv_sec) * 1000.0 +
                       (end.tv_usec - start.tv_usec) / 1000.0
                << std::endl;
        }
    }
    out.close();

}
