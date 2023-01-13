#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

nx = 100
ny = 100
c1 = 0.33
c2 = 0.33
file = f"{c1}-{c2}.txt"
data = np.loadtxt(file).reshape(nx,ny,3)

plt.imshow(data)
plt.show()
