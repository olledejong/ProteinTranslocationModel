import matplotlib.pyplot as plt
import numpy as np
import os
from model import averages_file, ref_trace_file

plt.style.use('seaborn-v0_8-dark')

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
    fig.suptitle("Cytosolic and nuclear volumes over time")
    ax1.set_xlabel('Cell cycle progression')
    ax1.grid(False)
    ax1.set_ylabel("Cytosolic volume", color='orange')
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
    fig.suptitle("Cytosolic and nuclear protein abundances over time")
    ax1.set_xlabel('Cell cycle progression')
    ax1.grid(False)
    ax1.set_ylabel("Cytosolic protein abundance", color='orange')
    ax1.plot(tspan / 100, y1, color='orange')
    ax1.tick_params(axis='y', labelcolor='orange')
    ax2 = ax1.twinx()
    ax2.grid(False)
    ax2.set_ylabel("Nuclear protein abundance", color='darkred')
    ax2.plot(tspan / 100, y2, color='darkred')
    ax2.tick_params(axis='y', labelcolor='darkred')
    save_figure(f"./output/{averages_file}/abundances.png")


def plot_concentration_ratio(final_tspan, c_con, n_con, con_ratio, kIn_base, kOut_base, kin_mp, kout_mp, ref_trace):
    cell_cycle_prog = final_tspan / 100  # convert x axis to cell cycle progression

    # calculate curve similarity
    exp_data = np.zeros((100, 2))
    exp_data[:, 0] = cell_cycle_prog
    exp_data[:, 1] = con_ratio

    num_data = np.zeros((100, 2))
    num_data[:, 0] = cell_cycle_prog
    num_data[:, 1] = ref_trace

    # mean absolute error
    mae = np.mean(np.abs(exp_data - num_data))
    mae = round(mae, 4)

    # mean squared error
    mse = np.mean(np.square(exp_data - num_data))
    mse = round(mse, 4)

    # fit polynomial through the data to make it look neater
    polfit = np.polyfit(cell_cycle_prog, con_ratio, 10)
    poly_y = np.polyval(polfit, cell_cycle_prog)

    plt.plot(cell_cycle_prog, c_con, c='darkred', lw=2)
    plt.title("Cytosolic protein concentration")
    plt.xlabel("Cell cycle progression")
    plt.ylabel("Concentration")
    save_figure(f"./output/{averages_file}/cyt_concentration.png")

    plt.plot(cell_cycle_prog, n_con, c='darkred', lw=2)
    plt.title("Nuclear protein concentration")
    plt.xlabel("Cell cycle progression")
    plt.ylabel("Concentration")
    save_figure(f"./output/{averages_file}/nuc_concentration.png")

    fig, ax = plt.subplots()
    ax.plot(cell_cycle_prog, poly_y, c='darkred', lw=4, alpha=0.8, label="Model prediction")
    ax.plot(cell_cycle_prog, ref_trace, c='grey', lw=2, alpha=0.6, label=f"{ref_trace_file.split('.')[0]} reference")
    plt.title(
        f"Nuclear to cytosolic protein concentration ratio\nParams: kIn base"
        f": {round(kIn_base, 6)}, kOut base: {round(kOut_base, 6)}, kIn mp: {kin_mp}, kOut mp: {kout_mp}"
        f"\nMean absolute error: {mae}, Mean squared error: {mse}"
    )
    plt.xlabel("Cell cycle progression")
    plt.ylabel("Ratio")
    plt.legend()
    save_figure(f"./output/{averages_file}/nc_concentration_ratio.png")


def plot_multiple_cycles(final_tspan, cyt_ab_cycles, nuc_ab_cycles, num_cycles):
    t_axis = np.linspace(0, final_tspan[-1] * num_cycles, len(cyt_ab_cycles))
    fig, ax1 = plt.subplots()
    fig.suptitle(f"Cytosolic and nuclear protein abundances over time (multiple cycles)")
    ax1.set_xlabel('Cell cycle progression')
    ax1.grid(False)
    ax1.set_ylabel("Cytosolic protein abundance", color='orange')
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
    save_figure(f"./output/{averages_file}/rates.png")
