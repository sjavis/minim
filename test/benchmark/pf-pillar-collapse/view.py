#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb

nx = 160
ny = 200


def map_colors(data, colors):
    assert data.shape[-1] == len(colors)
    colors = np.array([to_rgb(c) for c in colors])**2 # Square to linearise RGB intensities
    return np.sqrt(1 - np.dot(data, 1-colors).clip(0,1)) # 1-color so that empty=white


def plot(filename):
    data = np.loadtxt(filename).reshape((nx,ny,3))
    data = np.swapaxes(data, 0, 1) # Swap x-y axes for visualisation
    solid = data[:,:,0]
    liquid = data[:,:,1]
    gas = data[:,:,2]
    plt.imshow(map_colors(data, ['k','b','w']), origin='lower')
    plt.contour(solid, levels=[0.5], colors='k')
    plt.contour(np.ma.masked_where(solid>0.5, liquid-gas), levels=[0], colors='b')
    plt.gca().set_aspect('equal')
    plt.tight_layout()
    plt.show()


plot("outputs/res-0.001/st-100/s-80/ca-30/p-0.000000.txt")
