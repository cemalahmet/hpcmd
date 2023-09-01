#ifndef __ATOMS_H
#define __ATOMS_H

#include "types.h"

class Atoms {
  public:
    VectorX_t masses; // not used for the specific simulations.
    VectorX_t pot_energy;
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;

    Atoms(const Positions_t &p)
        : masses{p.cols()},
          pot_energy(p.cols()),
          positions{p},
          velocities{3, p.cols()},
          forces{3, p.cols()}
        {
        pot_energy.setZero();
        velocities.setZero();
        forces.setZero();
    }

    Atoms(const Positions_t &p, const Velocities_t &v)
        : masses{p.cols()},
          pot_energy(p.cols()),
          positions{p},
          velocities{v},
          forces{3, p.cols()} {
        assert(p.cols() == v.cols());
        pot_energy.setZero();
        forces.setZero();
    }

    Atoms(const Positions_t &p, const Velocities_t &v, const Forces_t &f)
        : masses{p.cols()},
          pot_energy(p.cols()),
          positions{p},
          velocities{v},
          forces{f} {
        assert(p.cols() == v.cols());
        assert(p.cols() == f.cols());
        pot_energy.setZero();
    }

    Atoms(const size_t nb_atoms)
        : masses{nb_atoms},
          pot_energy(nb_atoms),
          positions{3, nb_atoms},
          velocities{3, nb_atoms},
          forces{3, nb_atoms} {
        positions.setZero();
        velocities.setZero();
        forces.setZero();
        pot_energy.setZero();
    }

    size_t nb_atoms() const {
        return positions.cols();
    }

    void resize(const long new_size) {
        positions.conservativeResize(3, new_size);
        velocities.conservativeResize(3, new_size);
        forces.conservativeResize(3, new_size);
        pot_energy.conservativeResize(new_size);
        masses.conservativeResize(new_size);
    }
};

#endif