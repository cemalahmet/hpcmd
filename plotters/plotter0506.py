import matplotlib.pyplot as plt
import numpy as np
import os

def plot_execution_times(filepath1, filepath2):
    x_axis = []
    times05 = []
    times06 = []

    with open(filepath1, "r") as f:
        i = 0
        for line in f:
            if i % 2 == 0:
                x_axis.append(float(line))
            else:
                times05.append(float(line))
            i += 1

    with open(filepath2, "r") as f:
        i = 0
        for line in f:
            if i % 2 == 1:
                times06.append(float(line))
            i += 1

    plt.plot(x_axis, times05, color='red')
    plt.plot(x_axis, times06, color='green')

    plt.grid()

    path = "plots/milestone5-6(2).png"
    plt.legend(["LJ with Direct Summation", "LJ with Cutoff"])
    plt.xlabel('Number of Atoms')
    plt.ylabel('Simulation Time (ms)')
    plt.title("Simulation Time vs. Size of the System for the Two LJ Implementations")
    plt.savefig(path, bbox_inches='tight')


plot_execution_times("../milestones/05/out.dat", "../milestones/06/out.dat")