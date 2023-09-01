#include "headers/atoms.h"
#include "headers/lj_cutoff.h"
#include "headers/verlet.h"
#include "headers/xyz.h"
#include "headers/ducastelle.h"
#include "headers/mpi_support.h"
#include "headers/domain.h"

#include <cmath>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    std::ofstream traj("traj.xyz");
    std::ofstream out("out.dat");

    Domain domain(MPI_COMM_WORLD, {30, 30, 30},
                  {1, 1, MPI::comm_size(MPI_COMM_WORLD)}, {0, 0, 1});

    constexpr double mass = 20405.736652;
    constexpr double cutoff = 10.0;
    //constexpr double kB = 8.617333262e-5;

    auto [names, positions]{read_xyz("cluster_923.xyz")};
    Atoms atoms{positions};

    domain.enable(atoms); // divides into subdomains

    double total_time = 20000;
    double timestep = 1;

    long nb_steps = floor(total_time / timestep);

    NeighborList neighbor_list;

    for (int i = 0; i < nb_steps; i++) {
        verlet_step1(atoms, mass, timestep);

        domain.exchange_atoms(atoms);
        domain.update_ghosts(atoms, 2 * cutoff);

        neighbor_list.update(atoms, cutoff);
        ducastelle(atoms, neighbor_list, cutoff);

        verlet_step2(atoms, mass, timestep);

        // now calculate the total potential energy only using non-ghost atoms
        double local_pot = 0.0, local_kin = 0.0;
        for (int k = 0; k < domain.nb_local(); k++) {
            local_pot += atoms.pot_energy(k);
            local_kin += (mass * 0.5) * atoms.velocities.col(k).square().sum();
        }

        double global_kin{MPI::allreduce(local_kin, MPI_SUM, MPI_COMM_WORLD)};
        double global_pot{MPI::allreduce(local_pot, MPI_SUM, MPI_COMM_WORLD)};

        domain.disable(atoms);
        if (domain.rank() == 0) {
            out << i << std::endl;
            out << global_kin << std::endl;
            out << global_pot << std::endl;
            if (i % 10 == 0) write_xyz(traj, atoms);
        }
        domain.enable(atoms);
    }
    out.close();
    traj.close();
    MPI_Finalize();
}

