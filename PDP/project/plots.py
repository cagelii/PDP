import matplotlib.pyplot as plt
import numpy as np
from math import sqrt

def analytical(x,y):
    return x*(1-x)*y*(1-y)

def main():
    filename = 'output.txt'
    with open(filename, 'r') as file:
        data = file.read().split()
    n = int(sqrt(len(data)))
    h = 1/(n+1)
    print(f'n = {n}')

    data = np.array(data, dtype=float).reshape((n, n))

    x = np.linspace(0, 1, n+2)[1:-1]  # Exclude the boundaries
    y = np.linspace(0, 1, n+2)[1:-1]
    x, y = np.meshgrid(x, y)
    dataAna = analytical(x,y)
    diff = data-dataAna
    print(f'Norm-diff: {np.linalg.norm(diff)}')
    
    # Create a 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot the surface
    ax.plot_surface(x, y, data, cmap='viridis')

    # Add labels
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # Save the plot
    fig.savefig('solution.png')

    ax = fig.add_subplot(111, projection='3d')

    # Plot the surface
    ax.plot_surface(x, y, dataAna, cmap='viridis')

    # Add labels
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # Save the plot
    fig.savefig('solutionAna.png')

    #Speed up
    p = np.array([1,2,4,8])
    timeStrong = np.array([1.488083, 0.839849, 0.351128, 0.174050])
    timeWeak = np.array([0.084654, 0.085389, 0.087523])

    fig, ax = plt.subplots()
    ax.plot(p, timeStrong[0]/timeStrong, label="Experimental")
    ax.plot(p, p, label="Ideal")
    ax.set_ylabel("Speed up")
    ax.set_xlabel("Processes")
    ax.legend()
    fig.savefig("speed.png")

    p = np.array([1,4,9])
    fig, ax = plt.subplots()
    ax.plot(p, timeWeak[0]/timeWeak, label="Experimental")
    ax.plot(p, np.ones(3), label="Ideal")
    ax.set_ylabel("Efficiency")
    ax.set_xlabel("Processes")
    ax.legend()
    fig.savefig("eff.png")



main()