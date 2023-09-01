## High Performance Computing: Molecular Dynamics in C++

To execute this project, follow the guide below:
```bash
# To use MPI on remote system, run the following
module load compiler/gnu devel/cmake mpi/openmpi

# Download the repository
git clone https://github.com/cemalahmet/hpcmd.git
cd hpcmd

# Run the script to generate the Mackay Icosahedra for milestone 7
sh ./mackay/create_objects.sh

# Create the build directory
mkdir build
cd build

# Configure & compile 
cmake -DCMAKE_BUILD_TYPE=Release -DUSE_MPI=ON ..
make

# To run executables with MPI (milestones 8 and 9):
cd milestones/09
mpirun -n 4 ./09

# To run executables without MPI:
cd milestones/07
./07

# To run tests:
make test
```

To run the plotters after generating output, in the _hpcmd_ directory simply execute:
```bash
cd plotters
python3 plotter04.py 
```


