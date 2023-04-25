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


cell_vols, nuc_vols, nuc_surface_areas = load_data()
# interpolated vol/area functions for time dependent retrieval of values
t_range = np.arange(len(nuc_surface_areas))
a_func = interp1d(t_range, nuc_surface_areas, kind='linear', bounds_error=False)
cv_func = interp1d(t_range, cell_vols, kind='linear', bounds_error=False)
nv_func = interp1d(t_range, nuc_vols, kind='linear', bounds_error=False)


def main():
    kd = 0.2  # degradation rate for protein
    kc = 0.6  # synthesis rate of protein in cytosol
    kin = 0.8  # rate of translocation into nucleus
    kout = 0.2  # rate of translocation out of nucleus

    # initial conditions
    cp0 = 91.08659014056978
    np0 = 28.248611459856743

    # solve the ode up until the event
    tspan_before = np.linspace(0.0, 91, 200)
    ini_cond = np.array([cp0, np0])
    sols_before = odeint(dptd, ini_cond, tspan_before, args=(kd, kc, kin, kout))

    # simulate event (reduce nuclear abundance)
    perc_to_remove = (np.amax(nuc_vols) - nuc_vols[-1]) / np.amax(nuc_vols)  # reduction proportional to loss in nuc vol
    sols_after_event = np.array([
        sols_before[:, 0][-1],  # nothing happens to the cellular abundance at nuclear split
        sols_before[:, 1][-1] - (sols_before[:, 1][-1] * perc_to_remove)  # remove the calc. percentage from the ab.
    ])

    # simulate part after the event
    tspan_after = np.linspace(91, 98, 20)
    sols_after = odeint(dptd, sols_after_event, tspan_after, args=(kd, kc, kin, kout))

    # final values (concatenated)
    final_tspan = np.concatenate((tspan_before, tspan_after))
    final_sols = np.concatenate((sols_before, sols_after))

    # plotting
    plotting.plot_abundances(final_tspan, final_sols[:, 0], final_sols[:, 1])
    plotting.plot_volume_ratio(t_range, nuc_vols, cell_vols)
    plotting.plot_abundance_ratio(final_tspan, final_sols)


if __name__ == '__main__':
    main()
    sys.exit(0)
