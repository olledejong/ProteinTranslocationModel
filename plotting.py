import os
import numpy as np
import matplotlib.pyplot as plt
from model import averages_file, ref_trace_file

plt.style.use('seaborn-v0_8-dark')


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


def plot_volumes(tspan, cyt_vols, nuc_vols, flag):
    """
    Function that can be used to plot the volumes that are used in the model at any time. For example; the raw volumes
    over the cell cycle, and the adjusted volumes (after manual tweaking) to fit the model better.

    :param tspan:
    :param cyt_vols:
    :param nuc_vols:
    :param flag:
    :return:
    """
    fig, ax1 = plt.subplots()
    fig.suptitle(f"Cytosolic and nuclear volumes over the cell cycle ({flag} manual alterations)", y=0.95)
    ax1.set_xlabel('Cell cycle progression')
    ax1.grid(False)
    ax1.set_ylabel("Cytosolic volume", color='darkcyan')
    ax1.plot(tspan / 100, cyt_vols, c='darkcyan', lw=3, alpha=0.8)
    ax1.tick_params(axis='y', labelcolor='darkcyan')
    ax2 = ax1.twinx()
    ax2.grid(False)
    ax2.set_ylabel("Nuclear volume", color='darkred')
    ax2.plot(tspan / 100, nuc_vols, c='darkred', lw=3, alpha=0.8)
    ax2.tick_params(axis='y', labelcolor='darkred')
    save_figure(f"./output/{averages_file}/combined_volumes_{flag}.png")


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


def plot_concentration_prediction(t_span, con_ratio, kIn_base, kOut_base):
    cell_cycle_prog = t_span / 100  # convert x axis to cell cycle progression

    con_ratio = con_ratio / np.average(con_ratio)

    # fit polynomial through the data to make it look neater
    polfit = np.polyfit(cell_cycle_prog, con_ratio, 10)
    poly_y = np.polyval(polfit, cell_cycle_prog)

    fig, ax = plt.subplots()
    ax.plot(cell_cycle_prog, con_ratio, c='grey', lw=2, alpha=0.6, label="Raw prediction")
    ax.plot(cell_cycle_prog, poly_y, c='darkred', lw=4, alpha=0.8, label="Polynomial prediction")
    plt.title(
        f"Sfp1 N/C concentration ratio"
        f"\nParams: kIn base: {round(kIn_base, 6)}, kOut base: {round(kOut_base, 6)}"
    )
    plt.xlabel("Cell cycle progression")
    plt.ylabel("N/C concentration ratio")
    plt.legend()
    save_figure(f"./output/{averages_file}/nc_concentration_ratio.png")


def plot_prediction_vs_reference(t_span, c_con, n_con, con_ratio, kIn_base, kOut_base, kin_mp, kout_mp, ref_trace, normalize):
    cell_cycle_prog = t_span / 100  # convert x axis to cell cycle progression

    if normalize:
        con_ratio = con_ratio / np.average(con_ratio)
        ref_trace = ref_trace / np.average(ref_trace)

    mae, mse = get_similarity_measure(cell_cycle_prog, con_ratio, ref_trace)

    # fit polynomial through the data to make it look neater
    polfit = np.polyfit(cell_cycle_prog, con_ratio, 10)
    poly_y = np.polyval(polfit, cell_cycle_prog)

    plot_separate_concentrations(c_con, cell_cycle_prog, n_con)

    fig, ax = plt.subplots()
    ax.plot(cell_cycle_prog, poly_y, c='darkred', lw=4, alpha=0.8, label="Model prediction")
    ax.plot(cell_cycle_prog, ref_trace, c='grey', lw=2, alpha=0.6, label=f"{ref_trace_file.split('.')[0]} reference")
    plt.title(
        f"Nuclear to cytosolic protein concentration ratio\nParams: kIn base"
        f": {round(kIn_base, 6)}, kOut base: {round(kOut_base, 6)}"
    )
    plt.xlabel("Cell cycle progression")
    plt.ylabel("Ratio")
    plt.legend()
    save_figure(f"./output/{averages_file}/nc_concentration_ratio.png")


def get_similarity_measure(cell_cycle_prog, con_ratio, ref_trace):
    """
    Calculates a measure for the similarity of the predicted concentration ratio and the reference trace

    :param cell_cycle_prog:
    :param con_ratio:
    :param ref_trace:
    :return:
    """
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
    return mae, mse


def plot_separate_concentrations(c_con, cell_cycle_prog, n_con):
    """
    Plots the nuclear and cytosolic protein concentrations in separate plots. Can be used in addition to plotting the
    concentration ratio.

    :param c_con:
    :param cell_cycle_prog:
    :param n_con:
    :return:
    """
    plt.plot(cell_cycle_prog, c_con, c='darkred', lw=2, alpha=0.8)
    plt.title("Cytosolic protein concentration")
    plt.xlabel("Cell cycle progression")
    plt.ylabel("Concentration")
    save_figure(f"./output/{averages_file}/cyt_concentration.png")
    plt.plot(cell_cycle_prog, n_con, c='darkred', lw=2, alpha=0.8)
    plt.title("Nuclear protein concentration")
    plt.xlabel("Cell cycle progression")
    plt.ylabel("Concentration")
    save_figure(f"./output/{averages_file}/nuc_concentration.png")


def plot_rates(t_values, k_out_values, k_in_values):
    """
    In order to intuitively assess what rate(s) (curves) underlay the observed concentration ratio, this function
    plots the rates in one figure.

    :param t_values:
    :param k_out_values:
    :param k_in_values:
    :return:
    """
    to_cc_progr = [x / 100 for x in t_values]
    plt.scatter(to_cc_progr, k_out_values, s=2, label="nucl. export rate")
    plt.scatter(to_cc_progr, k_in_values, s=2, label="nucl. import rate")
    plt.title("Nuclear import and export rates over the cell cycle")
    plt.xlabel("Cell cycle progression")
    plt.ylabel("Rate")
    plt.legend(loc='best')
    save_figure(f"./output/{averages_file}/rates.png")
