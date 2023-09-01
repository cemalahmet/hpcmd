import matplotlib.pyplot as plt
import numpy as np
import os

def plot_execution_times(filepath):
    temperature = []
    energy = []

    with open(filepath, "r") as f:
        i = 0
        for line in f:
            if i % 2 == 0:
                temperature.append(float(line))
            else:
                energy.append(float(line))
            i += 1

    plt.plot(temperature, energy, color='red')

    plt.grid()

    path = "plots/milestone7.png"
    plt.xlabel('Temperature (K)')
    plt.ylabel('Energy (eV)')
    plt.title("Energy vs. Temperature")
    plt.savefig(path, bbox_inches='tight')


def plot_properties_vs_cluster_size():
    sizes = [147, 309, 561, 923, 1415, 2057, 2869, 3871, 5083, 6525]

    melting_temp = [429.24, 598.713, 690.19, 879.791, 924.301, 933.87, 941.108, 949.265, 947.385, 952.021]
    heat_cap = [0.0337087, 0.103128, 0.146451, 0.405466, 0.913482, 1.05046, 1.17026, 1.43584, 1.4395, 1.92114]
    latency = [5, 40, 80, 220, 440, 540, 800, 1320, 1680, 2100]

    plt.plot(sizes, melting_temp, color='red')
    plt.grid()
    path = "plots/milestone7-melt.png"
    plt.xlabel('Cluster Size')
    plt.ylabel('Melting Temperature (K)')
    plt.title("Melting Temperature vs. Cluster Size")
    plt.savefig(path, bbox_inches='tight')

    plt.clf()

    plt.plot(sizes, heat_cap, color='red')
    plt.grid()
    path = "plots/milestone7-hcap.png"
    plt.xlabel('Cluster Size')
    plt.ylabel('Heat Capacity (eV/K)')
    plt.title("Heat Capacity vs. Cluster Size")
    plt.savefig(path, bbox_inches='tight')

    plt.clf()

    plt.plot(sizes, latency, color='red')
    plt.grid()
    path = "plots/milestone7-latency.png"
    plt.xlabel('Cluster Size')
    plt.ylabel('Latent Heat (eV)')
    plt.title("Latent Heat vs. Cluster Size")
    plt.savefig(path, bbox_inches='tight')

plot_execution_times("../milestones/07/out.dat")
plot_properties_vs_cluster_size()
