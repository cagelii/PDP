import numpy as np
import matplotlib.pyplot as plt

threads = np.array([1, 2, 4, 8])
timeStrong = np.array([0.044463, 0.023079, 0.011814, 0.006674])
timesWeak = np.array([0.005893, 0.005710, 0.005733, 0.005908])

fig1, ax1 = plt.subplots(1,2)
ax1[0].plot(threads, timeStrong)
ax1[0].set_ylabel("Time (s)")
ax1[0].set_xlabel("Processes")
ax1[0].set_title("Strong scalability")

ax1[1].plot(threads, timesWeak)
ax1[1].set_ylabel("Time (s)")
ax1[1].set_xlabel("Processes")
ax1[1].set_title("Weak scalability")
fig1.tight_layout()
#fig1.suptitle("Timings", fontsize=16)
fig1.savefig("times.png")

speedStrong = timeStrong[0]/timeStrong
figStrong, axStrong = plt.subplots()
axStrong.plot(threads, speedStrong, label = "Experimental")
axStrong.plot(threads, threads, label="Ideal speed up")
axStrong.set_ylabel("Speed up")
axStrong.set_xlabel("Processes")
axStrong.legend()
figStrong.savefig("speed.png")
