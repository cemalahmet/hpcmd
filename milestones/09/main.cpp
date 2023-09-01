#include "headers/atoms.h"
#include "headers/berendsen_thermostat.h"
#include "headers/lj_cutoff.h"
#include "headers/verlet.h"
#include "headers/xyz.h"
#include "headers/ducastelle.h"
#include "headers/mpi_support.h"
#include "headers/domain.h"

#include <iostream>

    double normSquared(Velocities_t &velocities, int k) {
    return velocities.col(k).square().sum();
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    std::ofstream traj("traj.xyz");
    std::ofstream out("out.dat");

    // FOR SMALL WHISKER
    double domain_size = 142.08603661 - 0.72124892; // (largest z - smallest z)
    Domain domain(MPI_COMM_WORLD, {50, 50, domain_size},
                  {1, 1, MPI::comm_size(MPI_COMM_WORLD)}, {0, 0, 1});

    // FOR LARGE WHISKER
    // double domain_size = 287.77831781 - 0.72124892;
    // Domain domain(MPI_COMM_WORLD, {70, 70, domain_size},
    //                  {1, 1, MPI::comm_size(MPI_COMM_WORLD)}, {0, 0, 1});

    constexpr double mass = 20405.736652;
    constexpr double cutoff = 10.0;
    constexpr double kB = 8.617333262e-5;

    constexpr double relaxation_time = 250;
    constexpr double goal_temp = 0;

    auto [names, positions]{read_xyz("whisker_small.xyz")};
    Atoms atoms{positions};

    domain.enable(atoms);

    double timestep = 1;

    long scale_start = 1000;
    long scale_interval = 1000;
    double scale_amount = 0.5;

    // FOR SMALL WHISKER
    double area = 1600;

    // scale rate (unit = 1/fs)
    double scale_rate = scale_amount / (scale_interval * domain_size);
    std::cout << scale_rate << std::endl;

    double total_scale_amount = 20.0;
    double total_scaled = 0;

    NeighborList neighbor_list;

    int i = 0;
    while (total_scaled < total_scale_amount) {
        // the loop goes on until the nanowire is strained to the target value
        if (i >= scale_start && i % scale_interval == 0) {
            total_scaled += scale_amount;
            domain.scale(atoms, Vector3_t(50, 50, domain_size + total_scaled));
        }
        verlet_step1(atoms, mass, timestep);

        domain.exchange_atoms(atoms);
        domain.update_ghosts(atoms, 2 * cutoff);

        neighbor_list.update(atoms, cutoff);
        ducastelle(atoms, neighbor_list, cutoff);

        verlet_step2(atoms, mass, timestep);

        berendsen_thermostat(atoms, mass, goal_temp, timestep, relaxation_time, kB);

        // compute the forces on ghost atoms outside the boundaries
        double local_force1 = 0.0, local_force2 = 0.0;

        for (int j = domain.nb_local(); j < atoms.nb_atoms(); ++j) {
            // Ghost atoms to the bottom of the whisker
            if (atoms.positions(2, j) < 0.72124892)
                local_force1 += atoms.forces(2, j);

            // Ghost atoms to the top of the whisker
            if (atoms.positions(2, j) > domain_size + total_scaled)
                local_force2 += atoms.forces(2, j);
        }

        // Sum over all processes
        double global_force1{MPI::allreduce(local_force1, MPI_SUM, MPI_COMM_WORLD)};
        double global_force2{MPI::allreduce(local_force2, MPI_SUM, MPI_COMM_WORLD)};

        domain.disable(atoms);
        if (domain.rank() == 0) {
            if (i >= scale_start && i % scale_interval == 0) {
                double strain = total_scaled / domain_size;
                double stress = (global_force1 - global_force2) / area;
                out << strain << std::endl;
                out << stress << std::endl;
            }
            if (i % 10 == 0) write_xyz(traj, atoms);

            // To keep track of which iteration we are in
            if (!(i & (i-1))) std::cout << i << std::endl;
        }
        domain.enable(atoms);
        i++;
    }
    out.close();
    traj.close();
    MPI_Finalize();
}

