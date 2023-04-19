import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.interpolate import CubicSpline, interp1d, UnivariateSpline

averages = pd.read_excel("C:/Users/Olle de Jong/Documents/MSc Biology/MSB Research/Code and data/Data/Output/excel/final_averages_for_model.xlsx")

surface_areas = averages.cell_surface_areas.values
nuc_vols = averages.nuclear_volumes.values
cell_vols = averages.cell_volumes.values

timepoints = np.arange(len(surface_areas))
a_func = interp1d(timepoints, surface_areas, kind='cubic', bounds_error=False)
cv_func = interp1d(timepoints, cell_vols, kind='cubic', bounds_error=False)
nv_func = interp1d(timepoints, nuc_vols, kind='cubic', bounds_error=False)


def plot_abundances(tspan, y1, y2):
    fig, ax1 = plt.subplots()
    fig.suptitle("Cytosolic and nuclear protein abundances over time")
    ax1.set_xlabel('Time')
    ax1.grid(False)
    ax1.set_ylabel("Cytosolic protein abundance", color='tab:red')
    ax1.plot(tspan, y1, color='tab:red')
    ax2 = ax1.twinx()
    ax2.grid(False)
    ax2.set_ylabel("Nuclear protein abundance", color='tab:blue')
    ax2.plot(tspan, y2, color='tab:blue')
    plt.show()


def dptd(y, t, kd, kc, kImp, kExp):
    c = y[0]  # cellular abundance of protein
    n = y[1]  # nuclear abundance of protein

    # system equations, describe the change in cellular and nuclear abundances
    dC_dt = kc * 5 + a_func(t) * (-kImp * c / (cv_func(t) - nv_func(t)) + kExp * a_func(t) / nv_func(t)) - kd * c
    dN_dt = a_func(t) * (kImp * c / (cv_func(t) - nv_func(t)) - kExp * n / nv_func(t)) - kd * n

    return [dC_dt, dN_dt]


def main():
    # parameters
    kd = 0.2  # degradation rate for protein
    kc = 0.6  # synthesis rate of protein in cytosol
    kin = 0.8  # rate of translocation into nucleus
    kout = 0.2  # rate of translocation out of nucleus

    # initial conditions
    y0c = 194.9264937
    y0n = 79.43575304
    t0 = 0

    # final x-range and step to integrate over.
    t_final = 99  # final t value
    deltax = 0.05  # t-step

    # lists to store the results in
    t = [t0]
    sols_c = [y0c]
    sols_n = [y0n]
    t2 = t0

    # manually integrate at each time step, and check for event sign changes at each step
    while t2 <= t_final:  # stop integrating at t_final
        t1 = t[-1]
        t2 = round(t1 + deltax, 2)

        prev_sols = np.array([sols_c[-1], sols_n[-1]])

        # Event at t = 92, the nucleus is midway its division, so this is the point where we remove an equal part of
        # the protein abundance as well.
        if t2 == 92.0:
            perc_to_remove = (np.amax(nuc_vols) - nuc_vols[-1]) / np.amax(nuc_vols)

            print("Nuclear division took place. Removing a part (based on volume loss) of the nuclear protein abundance")
            prev_sols = np.array(
                [
                    sols_c[-1],
                    sols_n[-1] * perc_to_remove
                ]
            )

        y_new = odeint(dptd, prev_sols, [t1, t2], args=(kd, kc, kin, kout))  # integrate from t1,sol1 to t2,sol2
        t += [t2]

        sols_c += [y_new[-1][0]]
        sols_n += [y_new[-1][1]]

    plot_abundances(t, sols_c, sols_n)


if __name__ == '__main__':
    main()
    sys.exit(0)
