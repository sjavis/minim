import sys
import numpy as np

plot = True if (len(sys.argv)!=2) else bool(int(sys.argv[1]))


data = np.loadtxt("scaling.txt", skiprows=1)
n = data[:,0]
time = data[:,1]
mem = data[:,2]


if (plot):
    import matplotlib.pyplot as plt
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(n, time, 'ks-')
    ax2.plot(n, mem, 'rs-')
    ax1.set_xlabel("Num procs")
    ax1.set_ylabel("Time / s")
    ax2.set_ylabel("Memory / kB", color='r')
    plt.tight_layout()
    plt.show()
