import matplotlib.pyplot as plt


def plot_abundances(tspan, y1, y2):
    fig, ax1 = plt.subplots()
    fig.suptitle("Cytosolic and nuclear protein abundances over time")
    ax1.set_xlabel('Time')
    ax1.grid(False)
    ax1.set_ylabel("Cytosolic protein abundance", color='orange')
    ax1.plot(tspan, y1, color='orange')
    ax1.tick_params(axis='y', labelcolor='orange')
    ax2 = ax1.twinx()
    ax2.grid(False)
    ax2.set_ylabel("Nuclear protein abundance", color='darkred')
    ax2.plot(tspan, y2, color='darkred')
    ax2.tick_params(axis='y', labelcolor='darkred')
    plt.show()
