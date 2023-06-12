import matplotlib.pyplot as plt
import numpy as np
import os
from model import averages_file

# plt.style.use('seaborn-v0_8-dark')

# according to the volume analysis script, the average duration of one cycle is 71.35 minutes.
# time_scalar = 0.7135  # scalar based on average duration of cycle to scale back to real minute axis


def save_figure(path, bbox_inches='tight', dpi=300):
    """
    Custom function that lets you save a pyplot figure and creates the directory where necessary
    """
    directory = os.path.split(path)[0]
    filename = os.path.split(path)[1]
    if directory == '':
        directory = '.'

    if not os.path.exists(directory):
        os.makedirs(directory)

    save_path = os.path.join(directory, filename)

    # Actually save the figure
    plt.savefig(save_path, bbox_inches=bbox_inches, dpi=dpi)
    plt.close()


def plot_volumes(tspan, cell_vols, nuc_vols):
    cyt_vols = [i - j for i, j in zip(cell_vols, nuc_vols)]

    fig, ax1 = plt.subplots()
    fig.suptitle(f"Cytoplasmic and nuclear volumes over time")
    ax1.set_xlabel('Cell cycle progression')
    ax1.grid(False)
    ax1.set_ylabel("Cytoplasmic volume", color='orange')
    ax1.plot(tspan / 100, cyt_vols, color='orange')
    ax1.tick_params(axis='y', labelcolor='orange')
    ax2 = ax1.twinx()
    ax2.grid(False)
    ax2.set_ylabel("Nuclear volume", color='darkblue')
    ax2.plot(tspan / 100, nuc_vols, color='darkblue')
    ax2.tick_params(axis='y', labelcolor='darkblue')
    save_figure(f"./output/{averages_file}/combined_volumes.png")


def plot_abundances(tspan, y1, y2):
    fig, ax1 = plt.subplots()
    fig.suptitle(f"Cytoplasmic and nuclear protein abundances over time")
    ax1.set_xlabel('Cell cycle progression')
    ax1.grid(False)
    ax1.set_ylabel("Cytoplasmic protein abundance", color='orange')
    ax1.plot(tspan / 100, y1, color='orange')
    ax1.tick_params(axis='y', labelcolor='orange')
    ax2 = ax1.twinx()
    ax2.grid(False)
    ax2.set_ylabel("Nuclear protein abundance", color='darkred')
    ax2.plot(tspan / 100, y2, color='darkred')
    ax2.tick_params(axis='y', labelcolor='darkred')
    save_figure(f"./output/{averages_file}/abundances.png")


def plot_concentration_ratio(final_tspan, one_cycle_cyt, one_cycle_nuc, cv_func, nv_func):
    # get cell / nuc volumes for all timepoints
    cell_vols = [cv_func(t).flatten()[0] for t in final_tspan]
    nuc_vols = [nv_func(t).flatten()[0] for t in final_tspan]
    cyt_vols = [i - j for i, j in zip(cell_vols, nuc_vols)]

    # with the abundances, the event is distinctly happening over one time-step, this is not the case for the volumes.
    # this is because the volumes are first manipulated to represent nuclear division at one time-step, but after that
    # they are interpolated using the interp1d method. This causes a more transient simulation of that division.
    # therefore we need to slightly alter these to prevent plotting artifacts.
    nuc_vols[181:184] = [nuc_vols[180]] * len(nuc_vols[181:184])  # TODO make this dynamic (non-hardcoded indexes)
    cyt_vols[-2] = cyt_vols[-3]

    # calculate the concentrations
    c_con = [i / j for i, j in zip(one_cycle_cyt.tolist(), cyt_vols)]
    n_con = [i / j for i, j in zip(one_cycle_nuc.tolist(), nuc_vols)]

    con_ratio = [i / j for i, j in zip(n_con, c_con)]

    plt.plot(final_tspan / 100, c_con, c='darkred', lw=2)
    plt.title(f"Cytoplasmic protein concentration")
    plt.xlabel("Cell cycle progression")
    plt.ylabel("Concentration")
    save_figure(f"./output/{averages_file}/cyt_concentration.png")

    plt.plot(final_tspan / 100, n_con, c='darkred', lw=2)
    plt.title(f"Nuclear protein concentration")
    plt.xlabel("Cell cycle progression")
    plt.ylabel("Concentration")
    save_figure(f"./output/{averages_file}/nuc_concentration.png")

    plt.plot(final_tspan / 100, con_ratio, c='darkred', lw=2)
    plt.title(f"Nuclear to cytoplasmic protein concentration ratio")
    plt.xlabel("Cell cycle progression")
    plt.ylabel("Ratio")
    save_figure(f"./output/{averages_file}/nc_concentration_ratio.png")


def plot_multiple_cycles(final_tspan, cyt_ab_cycles, nuc_ab_cycles, num_cycles):
    t_axis = np.linspace(0, final_tspan[-1] * num_cycles, len(cyt_ab_cycles))
    fig, ax1 = plt.subplots()
    fig.suptitle(f"Cytoplasmic and nuclear protein abundances over time (multiple cycles)")
    ax1.set_xlabel('Cell cycle progression')
    ax1.grid(False)
    ax1.set_ylabel("Cytoplasmic protein abundance", color='orange')
    ax1.plot(t_axis / 100, cyt_ab_cycles, color='orange')
    ax1.tick_params(axis='y', labelcolor='orange')
    ax2 = ax1.twinx()
    ax2.grid(False)
    ax2.set_ylabel("Nuclear protein abundance", color='darkred')
    ax2.plot(t_axis / 100, nuc_ab_cycles, color='darkred')
    ax2.tick_params(axis='y', labelcolor='darkred')
    save_figure(f"./output/{averages_file}/multiple_cycles.png")


def plot_rates(ts, kouts, kins):
    to_cc_progr = [x / 100 for x in ts]
    plt.scatter(to_cc_progr, kouts, s=2, label="nucl. export rate")
    plt.scatter(to_cc_progr, kins, s=2, label="nucl. import rate")
    plt.title("Nuclear import and export rates over the cell cycle")
    plt.xlabel("Cell cycle progression")
    plt.ylabel("Rate")
    plt.legend(loc='best')
    plt.show()
