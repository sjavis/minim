#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

nx = 120
ny = 60
file = f"output.txt"
data = np.loadtxt(file).reshape(nx,ny)[:,1:-1]

plt.contourf(data.T, cmap='Blues')
plt.colorbar()
plt.gca().set_aspect('equal')
plt.show()
