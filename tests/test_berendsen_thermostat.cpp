#include "headers/berendsen_thermostat.h"
#include <gtest/gtest.h>

TEST(BerendsenThermostatTest, Convergence) {
    constexpr int test_count = 10;
    for (int _i = 0; _i < test_count; _i++) {
        constexpr int n = 20;

        constexpr double mass = 1.0;
        constexpr double goal_temperature = 300;
        constexpr double timestep = 0.001;
        constexpr double relaxation_time = 5.0;
        constexpr int nb_timesteps = 100000;
        constexpr double kB = 1.0;

        Eigen::Array3Xd positions(3, n), velocities(3, n);

        positions.setRandom();
        velocities.setRandom();

        Atoms atoms{positions, velocities};

        double last_temp = temperature(atoms, mass, kB);

        double temps_total = last_temp;

        // Prove that it gets closer to the goal temperature(the average of
        // temperatures over 50 time-steps is closer to the goal temperature)
        for (int i = 1; i < nb_timesteps; i++) {
            berendsen_thermostat(atoms, mass, goal_temperature, timestep,
                                 relaxation_time, kB);
            temps_total += temperature(atoms, mass, kB);
            if (i && !(i % 50)) {
                double avg = (temps_total / 50);

                EXPECT_LE(abs(goal_temperature - avg),
                          abs(goal_temperature - last_temp));
                last_temp = avg;
                // std::cout << last_temp << std::endl;
                temps_total = 0;
            }
        }

        // Expect that it converges
        EXPECT_NEAR(last_temp, goal_temperature, 1e-5);
    }
}