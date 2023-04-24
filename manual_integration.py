import sys
import plotting
import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.interpolate import CubicSpline, interp1d, UnivariateSpline

peak_of_nuc_vol = 80  # the latest local high of the nuclear volume
nuc_div_tp = 91  # simulated point at which nuclear division takes place


#########################
### Data preparations ###
#########################


def load_data():
    """
    Loads data from the averages file, calls the alter_data function and returns the outcome.
    :return:
    """
    file = "./averages.xlsx"
    averages = pd.read_excel(file)

    return alter_data(averages)


def alter_data(averages):
    """
    Takes the by the script exported averages and alters these based on various assumptions:
      - Nuclear division takes place approximately 10 time-steps before the end (which represents whole-cell splitting)
      - Nuclear division duration is very short, and thus assumed to be instant (split off at one time-point)
      - Nuclear volume and surf area stays constant from local high up until division point
    :param averages:
    :return:
    """

    # nuclear volume over time (altered to simulate instant nuclear division)
    nv = averages.nuc_volumes.values
    nv[peak_of_nuc_vol:nuc_div_tp] = nv[peak_of_nuc_vol]
    nv[nuc_div_tp:] = [nv[-1]] * len(nv[nuc_div_tp:])

    # nuclear surface area over time (altered such that it simulates a proportional decrease compared to nuclear volume)
    nsa = averages.nuc_surface_areas.values
    nsa[peak_of_nuc_vol:nuc_div_tp] = nsa[peak_of_nuc_vol]
    nsa[nuc_div_tp:] = [nsa[-1]] * len(nsa[nuc_div_tp:])

    # whole-cell volumes
    cv = averages.cell_volumes.values

    return cv, nv, nsa


###############################
### The model & simulations ###
###############################


def dptd(y, t, kd, kc, kImp, kExp):
    c, n = y  # cellular and nuclear abundance of protein

    # system equations, describe the change in cellular and nuclear abundances
    dC_dt = kc * 50 + a_func(t) * (-kImp * c / (cv_func(t) - nv_func(t)) + kExp * a_func(t) / nv_func(t)) - kd * c
    dN_dt = a_func(t) * (kImp * c / (cv_func(t) - nv_func(t)) - kExp * n / nv_func(t)) - kd * n

    return [dC_dt, dN_dt]


cell_vols, nuc_vols, nuc_surface_areas = load_data()  # load and alter the data
# interpolated vol/area functions for time dependent retrieval of values
t_range = np.arange(len(nuc_surface_areas))
a_func = interp1d(t_range, nuc_surface_areas, kind='cubic', bounds_error=False)
cv_func = interp1d(t_range, cell_vols, kind='cubic', bounds_error=False)
nv_func = interp1d(t_range, nuc_vols, kind='cubic', bounds_error=False)


def main():
    kd = 0.2  # degradation rate for protein
    kc = 0.6  # synthesis rate of protein in cytosol
    kin = 0.8  # rate of translocation into nucleus
    kout = 0.2  # rate of translocation out of nucleus

    # initial conditions
    cp0 = 91.08659014056978
    np0 = 28.248611459856743
    t0 = 0

    # final x-range and step to integrate over.
    t_final = 100  # final t value
    deltax = 0.05  # t-step


    # lists to store the results in
    t = [t0]
    sols_c = [cp0]
    sols_n = [np0]
    t2 = t0

    # manually integrate at each time step, and check for event sign changes at each step
    while t2 <= t_final:  # stop integrating at t_final
        t1 = t[-1]
        t2 = round(t1 + deltax, 2)

        prev_sols = np.array([sols_c[-1], sols_n[-1]])

        # Event at t = 91, the nucleus is midway its division, so this is the point where we remove an equal part of
        # the protein abundance as well.
        if t2 == float(nuc_div_tp):
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
