import sys
import numpy as np
import matplotlib.pyplot as plt

plot = True if (len(sys.argv)==1) else bool(int(sys.argv[1]))

acceptable_error = 100

surfacetension = 100
resolution = 1e-3
spacing = 40
contactangles = np.array([30, 50, 70, 90, 110, 130, 150])


def get_critical_pressure(spacing, contactangle, surfacetension, resolution):
    filename = f'outputs/res-{resolution}/st-{surfacetension}/s-{spacing}/ca-{contactangle}/critical_pressure.txt'
    with open(filename) as f:
        line = f.readline().split()
        return float(line[0]), float(line[2])

def get_analytic_pressure(spacing, contactangle, surfacetension, resolution):
    return -2 * surfacetension * np.cos(contactangle*np.pi/180) / (spacing*resolution)


pressures = np.zeros(len(contactangles))
errors = np.zeros(len(contactangles))
for i, ca in enumerate(contactangles):
    pressures[i], errors[i] = get_critical_pressure(spacing, ca, surfacetension, resolution)
residuals = pressures - get_analytic_pressure(spacing, contactangles, surfacetension, resolution)


if (plot):
    cas_analytic = np.linspace(0, 180, 100)
    pressures_analytic = get_analytic_pressure(spacing, cas_analytic, surfacetension, resolution)

    fig, axs = plt.subplots(2, 1, figsize=(6,3), sharex=True)
    axs[0].errorbar(contactangles, pressures, errors, fmt='ks', capsize=2, mec='r', mfc='r')
    axs[0].plot(cas_analytic, pressures_analytic, 'k-')
    axs[1].errorbar(contactangles, residuals, errors, fmt='ks', capsize=2, mec='r', mfc='r')
    axs[1].axhspan(-acceptable_error, acceptable_error, color='green', alpha=0.2)
    axs[1].axhline(0, c='k', ls=':', lw=1)
    axs[0].set_ylabel("Critical Pressure")
    axs[1].set_ylabel("Residuals")
    axs[1].set_xlabel("Contact Angle")
    plt.tight_layout()
    plt.show()


if (np.any(np.abs(residuals) > acceptable_error)):
    sys.exit(1)
