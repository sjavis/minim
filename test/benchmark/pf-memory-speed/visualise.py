#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

nx = 200
ny = 200
file = f"output.txt"
data = np.loadtxt(file).reshape(nx,ny,3)

plt.imshow(data)
plt.show()
