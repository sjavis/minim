import sys
import numpy as np
from scipy.optimize import curve_fit

assert (len(sys.argv) >= 2)
widths = np.fromstring(sys.argv[1], sep=' ')
plot = True if (len(sys.argv)!=3) else bool(int(sys.argv[2]))


def interface(x, x0, width):
    return 0.5 + 0.5*np.tanh((x-x0)/(2*width))

def fit_interface(phi):
    x = np.arange(len(phi))
    [x0, width], _ = curve_fit(interface, x, phi, [len(phi)/2, 1])
    return x0, width


errors = np.zeros(len(widths))
for i_w, w0 in enumerate(widths):
    data = np.loadtxt(f"width-{w0:g}.txt")
    data = data.reshape((-1, 2))

    phi = data[:,0]
    x0, width = fit_interface(phi)
    errors[i_w] = np.abs(width - w0) / w0


if (plot):
    import matplotlib.pyplot as plt
    plt.plot(widths, errors, 'ks-')
    plt.loglog()
    plt.xlabel("Interface width")
    plt.ylabel("Error in measured width")
    plt.show()

if (errors[widths==1] > 0.025):
    sys.exit(1)
