#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

nc = 3


def load(file):
    data = np.loadtxt(file)
    nxy = int(np.sqrt(data.size/nc))
    return data.reshape((nxy,nxy,nc))


def plot(filename):
    data = load(filename)
    data = data.swapaxes(0,1)
    plt.imshow(data, origin='lower', cmap='Blues', vmin=0, vmax=1)
    solid = data[:,:,0]
    liquid = data[:,:,1]
    gas = data[:,:,2]
    plt.contour(solid, levels=[0.5], colors='k')
    plt.contour(np.ma.masked_where(solid>0.5, liquid-gas), levels=[0], colors='k')


plot('minimum.txt')

plt.tight_layout()
plt.show()
