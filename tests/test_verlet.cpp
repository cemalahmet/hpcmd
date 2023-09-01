#include "headers/verlet.h"
#include <gtest/gtest.h>

TEST(VerletTest, UnderConstantForce) {
    long nb_steps = 10000;
    double timestep = 0.01;
    double mass = 1.0;

    double total_time = timestep * nb_steps;

    int n = 10;
    Eigen::Array3Xd positions(3, n), velocities(3, n), forces(3, n);

    positions.setRandom();
    velocities.setRandom();
    forces.setRandom();

    Eigen::Array3Xd r_positions(3, n), r_velocities(3, n);
    for (int i = 0; i < n; i++) {
        auto a = forces.col(i);
        auto v = velocities.col(i);
        auto x = positions.col(i);

        for (int j = 0; j < 3; j++) {
            r_positions(j, i) = 0.5 * a(j) * total_time * total_time + v(j) * total_time + x(j);
            r_velocities(j, i) = a(j) * total_time + v(j);
        }

    }

    for (int i = 0; i < nb_steps; ++i) {
        verlet_step1(positions, velocities, forces, mass, timestep);
        verlet_step2(velocities, forces, mass, timestep);
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 3; j++) {
            EXPECT_NEAR(velocities(j,i), r_velocities(j,i), 1e-5);
            EXPECT_NEAR(positions(j,i), r_positions(j,i), 1e-5);
        }
    }
}