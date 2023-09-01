#include "headers/hello.h"
#include <Eigen/Dense>
#include <filesystem>
#include <iostream>

#include <mpi.h>

int main(int argc, char *argv[]) {
    int rank = 0, size = 1;

    // Below is some MPI code, try compiling with `cmake -DUSE_MPI=ON ..`
    MPI_Init(&argc, &argv);

    // Retrieve process infos
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::cout << "Hello I am rank " << rank << " of " << size << "\n";

    if (rank == 0)
      hello_eigen();

    auto input_path = "./simulation_test_input.txt";

    if (not std::filesystem::exists(input_path))
      std::cerr << "warning: could not find input file " << input_path << "\n";

    MPI_Finalize();

    return 0;
}
