#ifndef __BERENDSEN_THERMOSTAT_H
#define __BERENDSEN_THERMOSTAT_H

#include "atoms.h"

double temperature(Atoms &atoms, double mass, double kB);
double temperature(double Ek, size_t n, double kB);

void berendsen_thermostat(Atoms &atoms, double mass, double goal_temperature,
                          double timestep, double relaxation_time, double kB);

#endif // __BERENDSEN_THERMOSTAT_H
