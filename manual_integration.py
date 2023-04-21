import sys
import plotting
import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.interpolate import CubicSpline, interp1d, UnivariateSpline

file = "C:/Users/Olle de Jong/Documents/MSc Biology/MSB Research/Code and data/Data/Output/excel/final_averages_for_model.xlsx"
averages = pd.read_excel(file)

peak_of_nuc_vol = 80

# nuclear volume over time (altered to simulate instant nuclear division)
nuc_vols = averages.nuc_volumes.values
nuc_vols[peak_of_nuc_vol:] = [nuc_vols[-1]] * len(nuc_vols[peak_of_nuc_vol:])

# nuclear surface area over time (altered such that it simulates a proportional decrease compared to nuclear volume)
nuc_surface_areas = averages.nuc_surface_areas.values
nuc_surface_areas[peak_of_nuc_vol:] = [nuc_surface_areas[-1]] * len(nuc_surface_areas[peak_of_nuc_vol:])

# whole-cell volumes
cell_vols = averages.cell_volumes.values

# interpolated vol/area functions for time dependent retrieval of values
timepoints = np.arange(len(nuc_surface_areas))
a_func = interp1d(timepoints, nuc_surface_areas, kind='cubic', bounds_error=False)
cv_func = interp1d(timepoints, cell_vols, kind='cubic', bounds_error=False)
nv_func = interp1d(timepoints, nuc_vols, kind='cubic', bounds_error=False)


def dptd(y, t, kd, kc, kImp, kExp):
    c, n = y  # cellular and nuclear abundance of protein

    # system equations, describe the change in cellular and nuclear abundances
    dC_dt = kc * 50 + a_func(t) * (-kImp * c / (cv_func(t) - nv_func(t)) + kExp * a_func(t) / nv_func(t)) - kd * c
    dN_dt = a_func(t) * (kImp * c / (cv_func(t) - nv_func(t)) - kExp * n / nv_func(t)) - kd * n

    return [dC_dt, dN_dt]


def main():
    # parameters
    kd = 0.2  # degradation rate for protein
    kc = 0.6  # synthesis rate of protein in cytosol
    kin = 0.8  # rate of translocation into nucleus
    kout = 0.2  # rate of translocation out of nucleus

    # initial conditions
    y0c = 91.4433196240664
    y0n = 28.571177299490948
    t0 = 0

    # final x-range and step to integrate over.
    t_final = 99.95  # final t value
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
        if t2 == float(peak_of_nuc_vol):
            print("Nuclear division took place. Removing a part (based on volume loss) of the nuclear protein abundance")

            perc_to_remove = (np.amax(nuc_vols) - nuc_vols[-1]) / np.amax(nuc_vols)
            prev_sols = np.array([sols_c[-1], sols_n[-1] - (sols_n[-1] * perc_to_remove)])

        y_new = odeint(dptd, prev_sols, [t1, t2], args=(kd, kc, kin, kout))  # integrate from t1,sol1 to t2,sol2

        t += [t2]

        sols_c += [y_new[-1][0]]
        sols_n += [y_new[-1][1]]

    sols_c = [x for x in sols_c if str(x) != 'nan']
    sols_n = [x for x in sols_n if str(x) != 'nan']

    plotting.plot_abundances(t[:len(sols_c)], sols_c, sols_n)


if __name__ == '__main__':
    main()
    sys.exit(0)
