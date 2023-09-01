#include "headers/atoms.h"
#include "headers/lj_direct_summation.h"
#include "headers/lj_cutoff.h"
#include "headers/verlet.h"
#include "headers/xyz.h"
#include "headers/berendsen_thermostat.h"

#include <iostream>
#include <cmath>
#include <sys/time.h>

int main(int argc, char *argv[]) {

    std::ofstream traj("traj.xyz");
    std::ofstream out("out.dat");

    constexpr double sigma = 1.0;
    constexpr double epsilon = 1.0;
    constexpr double mass = 1.0;
    constexpr double cutoff = 4.0;

    double lattice_constant{1.4 * sigma};

    double total_time = 100 * sigma * sqrt(mass / epsilon);
    double timestep = 0.002 * sigma * sqrt(mass / epsilon);

    double goal_temperature = 0.35;
    double kB = 1.0;

    long nb_steps = floor(total_time / timestep);

    struct timeval start, end;

    // Test for performance:
    // the structures tested are 2*2*2, 2*2*3, 2*3*3, 3*3*3, 3*3*4, ... , 6*7*7
    for (int N = 2; N <= 6; N++) {
        for (int M = 0; M < 3; M++) {
            int N1 = N;
            int N2 = (N + (M > 0));
            int N3 = (N + (M > 1));

            Atoms atoms{size_t(N1 * N2 * N3)};

            // Cubic lattice, small random dislocations
            for (int x{0}, i{0}; x < N1; ++x) {
                for (int y{0}; y < N2; ++y) {
                    for (int z{0}; z < N3; ++z, ++i) {
                        atoms.positions(0, i) += x * lattice_constant;
                        atoms.positions(1, i) += y * lattice_constant;
                        atoms.positions(2, i) += z * lattice_constant;
                    }
                }
            }

            NeighborList neighbor_list;
            gettimeofday(&start, NULL);
            for (int i = 0; i < nb_steps; i++) {
                verlet_step1(atoms, mass, timestep);

                neighbor_list.update(atoms, cutoff);
                lj_cutoff(atoms, neighbor_list, epsilon, sigma, cutoff);

                verlet_step2(atoms, mass, timestep);


                berendsen_thermostat(atoms, mass, goal_temperature, timestep,
                                         100 * timestep, kB);
            }
            gettimeofday(&end, NULL);

            out << N1 * N2 * N3 << std::endl;

            // output in miliseconds
            out << (end.tv_sec - start.tv_sec) * 1000.0 +
                       (end.tv_usec - start.tv_usec) / 1000.0
                << std::endl;
        }
    }
    out.close();
    traj.close();
}
