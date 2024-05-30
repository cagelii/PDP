import numpy as np
import matplotlib.pyplot as plt

timesstrong = np.array([26.04, 13.5, 7.09, 3.92])
timesweak = np.array([26.04, 28.04, 30.96, 36.38])
processes = np.array([1, 2, 4, 8])

speedstrong = timesstrong[0]/timesstrong

fig, ax  = plt.subplots()
fig.suptitle("Speed up for strong scalability")
ax.plot(processes, speedstrong, label = "Experimental")
ax.plot(processes, processes, label = "Ideal")
ax.set_ylabel("Speed up")
ax.set_xlabel("Processes")
ax.legend()
fig.savefig("plotstrong.png")

fig, ax  = plt.subplots()
fig.suptitle("Efficiency for weak scalability")
ax.plot(processes, timesweak[0]/timesweak, label = "Experimental")
ax.set_ylabel("Efficiency")
ax.set_xlabel("Processes")
fig.savefig("plotweak.png")
