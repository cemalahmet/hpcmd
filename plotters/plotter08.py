import matplotlib.pyplot as plt
import numpy as np
import os

def plot_energies(filepath):
    with open(filepath, "r") as f:
        i = 0
        kin = []
        pot = []
        tot = []
        for line in f:
            if i % 3 == 0:
                i += 1
                continue
            elif i % 3 == 1:
                i += 1
                kin.append(float(line))
            else:
                i += 1
                pot.append(float(line))
                tot.append(kin[-1] + pot[-1])

        timestep = 1
        x_axis = np.linspace(timestep, len(kin) * timestep, len(kin))

        plt.plot(x_axis, kin, color='red')
        plt.plot(x_axis, pot, color='green')
        plt.plot(x_axis, tot, color='blue')
        plt.grid()

        path = "plots/milestone8.png"
        plt.legend(["Kinetic Energy", "Potential Energy", "Total Energy"])
        plt.xlabel('Time (fs)')
        plt.ylabel('Energy (eV)')
        plt.title("Conservation of Energy Simulation on Gold Cluster")
        plt.savefig(path, bbox_inches='tight')

        plt.clf()

        plt.plot(x_axis[::100], tot[::100], color='blue')
        plt.ylim([-3316, -3318])
        plt.xlim([-1000, 21000])
        plt.grid()
        path = "plots/milestone8-only-total.png"
        plt.xlabel('Time (fs)')
        plt.ylabel('Total Energy (eV)')
        plt.title("Total Energy vs. Time of a Gold Cluster")
        plt.savefig(path, bbox_inches='tight')

plot_energies("../milestones/08/out.dat")
