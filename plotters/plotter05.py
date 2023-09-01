import matplotlib.pyplot as plt
import numpy as np

def plot_temp(filepath1, filepath2, filepath3):
    timestep = 0.002
    temp1 = []
    temp2 = []
    temp3 = []
    plot_period = 150

    with open(filepath1, "r") as f:
        i = 0
        for line in f:
            if i % plot_period == 0:
                temp1.append(float(line))
            i += 1


    with open(filepath2, "r") as f:
        i = 0
        for line in f:
            if i % plot_period == 0:
                temp2.append(float(line))
            i += 1


    with open(filepath3, "r") as f:
        i = 0
        for line in f:
            if i % plot_period == 0:
                temp3.append(float(line))
            i += 1

    x_axis = np.linspace(timestep, len(temp1) * timestep * plot_period, len(temp1))

    plt.plot(x_axis, temp1, color='red')
    plt.grid()

    path = "plots/milestone5_1.png"
    plt.xlabel('Time (LJ reduced units)')
    plt.ylim([-0.05, 1.05])
    plt.ylabel('Temperature (LJ reduced units)')
    plt.title("Temperature vs. Time")
    plt.savefig(path, bbox_inches='tight')

    plt.clf()
    plt.plot(x_axis, temp2, color='green')
    plt.grid()

    path = "plots/milestone5_2.png"
    plt.xlabel('Time (LJ reduced units)')
    plt.ylim([-0.05, 1.05])
    plt.ylabel('Temperature (LJ reduced units)')
    plt.title("Temperature vs. Time")
    plt.savefig(path, bbox_inches='tight')

    plt.clf()
    plt.plot(x_axis, temp3, color='blue')
    plt.grid()

    path = "plots/milestone5_3.png"
    plt.xlabel('Time (LJ reduced units)')
    plt.ylim([-0.05, 1.05])
    plt.ylabel('Temperature (LJ reduced units)')
    plt.title("Temperature vs. Time")
    plt.savefig(path, bbox_inches='tight')


plot_temp("../milestones/05/out1.dat","../milestones/05/out2.dat", "../milestones/05/out3.dat")
