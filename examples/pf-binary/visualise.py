#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

nx = 100
ny = 100
file = f"output.txt"
data = np.loadtxt(file).reshape(nx,ny)

plt.contourf(data.T, cmap='Blues')
plt.colorbar()
plt.show()
