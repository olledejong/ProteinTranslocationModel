import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.interpolate import CubicSpline, interp1d, UnivariateSpline

averages = pd.read_excel("C:/Users/Olle de Jong/Documents/MSc Biology/MSB Research/Code and data/Data/Output/excel/averages.xlsx")

surface_areas = averages.cell_surface_areas.values
nuc_vols = averages.cell_volumes.values
cell_vols = averages.nuclear_volumes.values

timepoints = np.arange(len(surface_areas))
a_func = interp1d(timepoints, surface_areas, kind='cubic', bounds_error=False)
cv_func = interp1d(timepoints, nuc_vols, kind='cubic', bounds_error=False)
nv_func = interp1d(timepoints, cell_vols, kind='cubic', bounds_error=False)


# TODO implement what happens to abundances at nuclear division. At t = 92, the nuclear division is half way done
# TODO at that moment the abundance in the nucleus should decrease. This needs to be ~ 26%
def dptd(y, t, kd, kc, kImp, kExp):
    c = y[0]  # cellular abundance of protein
    n = y[1]  # nuclear abundance of protein

    # system equations, describe the change in cellular and nuclear abundances
    dC_dt = kc * 5 + a_func(t) * (-kImp * c / (cv_func(t) - nv_func(t)) + kExp * a_func(t) / nv_func(t)) - kd * c
    dN_dt = a_func(t) * (kImp * c / (cv_func(t) - nv_func(t)) - kExp * n / nv_func(t)) - kd * n

    return [dC_dt, dN_dt]


def main():
    kd = 0.2  # degradation rate for protein
    kc = 0.6  # synthesis rate of protein in cytosol
    kin = 0.8  # rate of translocation into nucleus
    kout = 0.2  # rate of translocation out of nucleus

    tspan = np.linspace(0.0, 99, 200)  # timepoints for which we want the output
    ini_cond = np.array([194.9264937, 79.43575304])  # initial abundances
    print(tspan)

    # Solving the ODE
    yt = odeint(
        func=dptd,
        y0=ini_cond,
        t=tspan,
        args=(kd, kc, kin, kout)
    )
    plot_abundances(tspan, yt)


def plot_abundances(tspan, yt):
    fig, ax1 = plt.subplots()
    fig.suptitle("Cytosolic and nuclear protein abundances over time")
    ax1.set_xlabel('Time')
    ax1.grid(False)
    ax1.set_ylabel("Cytosolic protein abundance", color='tab:red')
    ax1.plot(tspan, yt[:, 0], color='tab:red')
    ax1.tick_params(axis='y', labelcolor='tab:red')
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax2.grid(False)
    ax2.set_ylabel("Nuclear protein abundance", color='tab:blue')
    ax2.plot(tspan, yt[:, 1], color='tab:blue')
    ax2.tick_params(axis='y', labelcolor='tab:blue')
    plt.show()


if __name__ == '__main__':
    main()
    sys.exit(0)
