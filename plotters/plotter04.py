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
                pot.append(float(line))
            else:
                i += 1
                kin.append(float(line))
                tot.append(kin[-1] + pot[-1])

        timestep = 0.001
        x_axis = np.linspace(timestep, len(kin) * timestep, len(kin))

        plt.plot(x_axis, kin, color='red')
        plt.plot(x_axis, pot, color='green')
        plt.plot(x_axis, tot, color='blue')
        plt.grid()

        path = "plots/milestone4.png"
        plt.legend(["Kinetic Energy", "Potential Energy", "Total Energy"])
        plt.xlabel('Time (LJ reduced units)')
        plt.ylabel('Energy (LJ reduced units)')
        plt.title("Conservation of Energy Simulation")
        plt.savefig(path, bbox_inches='tight')


def plot_energy_vs_timestep(filepath):
    with open(filepath, "r") as f:
        en = [[]]
        for line in f:
            if line[0] == '#':
                en.append([])
                continue
            else:
                en[-1].append(float(line))

        ts = [0.001, 0.002, 0.005, 0.01, 0.02]
        x_axis = np.linspace(0.02,100,5000)

        plt.plot(x_axis,en[0][::20],  color='red')
        plt.plot(x_axis,en[1][::10],  color='orange')
        plt.plot(x_axis,en[2][::4],  color='green')
        plt.plot(x_axis,en[3][::2],  color='blue')
        plt.plot(x_axis,en[4],       color='purple')

        plt.grid()

        path = "plots/milestone4_2.png"
        plt.legend(["0.001", "0.002", "0.005", "0.01", "0.02"])
        plt.xlabel('Time (LJ reduced units)')
        plt.ylabel('Total Energy (LJ reduced units)')
        plt.title("Total Energy Over Time for Different Timesteps")
        plt.savefig(path, bbox_inches='tight')

        for energy_list in en:
            if energy_list:
                print(np.std(energy_list))

        for energy_list in en:
            if energy_list:
                print(np.average(energy_list))

plot_energies("../milestones/04/out.dat")
plt.clf()
plot_energy_vs_timestep("../milestones/04/out2.dat")