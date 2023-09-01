import matplotlib.pyplot as plt
import numpy as np
import os

def plot_stress_strain(filepath):
    with open(filepath, "r") as f:
        i = 0

        # add the point (0,0) by default
        strain = [0]
        stress = [0]

        for line in f:
            if i % 2 == 0:
                strain.append(float(line))
                i += 1
            elif i % 2 == 1:
                stress.append(float(line))
                i += 1

        plt.plot(strain, stress, color='red')
        plt.grid()

        path = "plots/milestone9_6.png"
        plt.xlabel('Strain')
        plt.ylabel('Stress (eV/Angstrom^3)')
        plt.ylim([-0.001, 0.016])
        plt.title("Stress vs. Strain")
        plt.savefig(path, bbox_inches='tight')

plot_stress_strain("../milestones/09/out6.dat")
