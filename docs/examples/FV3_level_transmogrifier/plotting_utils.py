import math
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate as sci_interp

DASH_PATTERNS = [
    [5, 0, 5, 0],
    [5, 2, 5, 2],
    [5, 2, 2, 2],
    [2, 2, 2, 2],
    [10, 3, 10, 3],
]
LINE_COLORS = ["black", "red", "darkgreen", "orange", "gray"]
MARKERSIZE = 5
ZLIMIT_CHOICES = ["PBL"]


def get_best_subplot_layout(num_plots):
    """
    Calculates the optimal subplot layout (rows, cols) for a given number of plots.
    """

    # If only one plot, 1x1 is the obvious choice
    if num_plots == 1:
        return 1, 1

    # Find the two factors of num_plots that are closest together
    sqrt_plots = math.sqrt(num_plots)
    rows = int(sqrt_plots)
    cols = math.ceil(num_plots / rows)

    # If the product of rows and cols is less than num_plots,
    # increment one of them to accommodate all plots
    while rows * cols < num_plots:
        if rows < cols:
            rows += 1
        else:
            cols += 1

    return rows, cols


def get_subplot_indices(plt_count, nrows, ncols):
    iy = math.floor(plt_count / ncols)
    ix = plt_count - iy * ncols
    return iy, ix


def fv3_level_make_plots(
    ab_dict, zss, do_ylog, zcoord_type, zlimit_type=None, save=False
):
    if not zlimit_type is None:
        assert zlimit_type in ZLIMIT_CHOICES, "Not a valid choice for zlimit_type"
    # Load standard atmosphere
    # Taken from /home/lmh/research/xic_levels

    DPI_FOR_PLOT = 300

    data = np.loadtxt("std_atmos_1976.txt")
    zstd = data[:, 0]
    Tstd = data[:, 1]
    pstd = data[:, 2]

    # Prepare interpolation between z and p
    Iz2p = sci_interp.interp1d(zstd, pstd, kind="linear")
    Ip2z = sci_interp.interp1d(pstd, zstd, kind="linear")

    num_plots = 7
    nrows, ncols = get_best_subplot_layout(num_plots)
    fig, ax = plt.subplots(nrows, ncols, figsize=(24, 11), sharey=True)

    assert len(DASH_PATTERNS) >= len(ab_dict), "Not enough dash patterns defined."

    for count, config in enumerate(ab_dict):
        line_color = LINE_COLORS[count]
        ak = ab_dict[config]["ak"]
        bk = ab_dict[config]["bk"]

        for zs_count, zs in enumerate(zss):
            dash_pattern = DASH_PATTERNS[zs_count]
            # Compute p and auxiliary quantities
            ps = Iz2p(zs)
            pe = ak + bk * ps
            delp = np.diff(pe)
            logpe = np.log(pe)
            dlogp = np.diff(logpe)
            pm = delp / dlogp

            ze = Ip2z(pe)
            dz = -np.diff(ze)
            zm = 0.5 * (ze[:-1] + ze[1:])

            npz = pm.size
            ke = np.arange(npz + 1) - 0.5
            km = np.arange(npz)

            zcoords_m = {"p": pm, "z": zm, "k": km}
            zcoords_e = {"p": pe, "z": ze, "k": ke}

            zcoord_labels = {
                "p": "Pressure [Pa]",
                "z": "Height ASL [m]",
                "k": "Coordinate Layer",
            }
            zlimits = {
                "p": {"PBL": (100000, 85000)},
                "z": {"PBL": (0, 4000)},
                "k": {"PBL": (30, 50)},
            }

            zcoord_m = zcoords_m[zcoord_type]
            zcoord_e = zcoords_e[zcoord_type]

            plt_count = 0
            label = config + " " + f"{zs/1000.:.0f}" + " km"
            ax[get_subplot_indices(plt_count, nrows, ncols)].plot(
                delp,
                zcoord_m,
                dashes=dash_pattern,
                markersize=MARKERSIZE,
                label=label,
                color=line_color,
                marker="o",
            )
            ax[get_subplot_indices(plt_count, nrows, ncols)].set_title("delp [Pa]")
            if zcoord_type != "k":
                dxax = np.diff(
                    ax[get_subplot_indices(plt_count, nrows, ncols)].get_xlim()
                )[0]
                ax[get_subplot_indices(plt_count, nrows, ncols)].text(
                    0 + dxax * (0.05 + count * 0.2),
                    zcoord_m[-1],
                    "p_s = %d hPa\nz_s = %d m" % (ps / 100.0, zs),
                    horizontalalignment="left",
                    verticalalignment="center",
                    color=line_color,
                )

            plt_count += 1
            ax[get_subplot_indices(plt_count, nrows, ncols)].plot(
                delp[:-1] / delp[1:],
                zcoord_e[1:-1],
                dashes=dash_pattern,
                markersize=MARKERSIZE,
                color=line_color,
                marker="o",
            )
            ax[get_subplot_indices(plt_count, nrows, ncols)].set_title("ratio delp")

            plt_count += 1
            ax[get_subplot_indices(plt_count, nrows, ncols)].plot(
                dlogp,
                zcoord_m,
                dashes=dash_pattern,
                markersize=MARKERSIZE,
                color=line_color,
                marker="o",
            )
            ax[get_subplot_indices(plt_count, nrows, ncols)].set_title("dlogp [Pa]")

            plt_count += 1
            ax[get_subplot_indices(plt_count, nrows, ncols)].plot(
                dlogp[:-1] / dlogp[1:],
                zcoord_e[1:-1],
                dashes=dash_pattern,
                markersize=MARKERSIZE,
                color=line_color,
                marker="o",
            )
            ax[get_subplot_indices(plt_count, nrows, ncols)].set_title("ratio dlogp")

            plt_count += 1
            ax[get_subplot_indices(plt_count, nrows, ncols)].plot(
                dz,
                zcoord_m,
                dashes=dash_pattern,
                markersize=MARKERSIZE,
                color=line_color,
                marker="o",
            )
            ax[get_subplot_indices(plt_count, nrows, ncols)].set_title("dz[m]")
            if zlimit_type == "PBL":
                ax[get_subplot_indices(plt_count, nrows, ncols)].set_xlim(0, 1000)

            plt_count += 1
            ax[get_subplot_indices(plt_count, nrows, ncols)].plot(
                dz[:-1] / dz[1:],
                zcoord_e[1:-1],
                dashes=dash_pattern,
                markersize=MARKERSIZE,
                color=line_color,
                marker="o",
            )
            ax[get_subplot_indices(plt_count, nrows, ncols)].set_title("ratio dz")

            plt_count += 1
            ax[get_subplot_indices(plt_count, nrows, ncols)].plot(
                np.linspace(npz, 1, npz),
                zcoord_m,
                dashes=dash_pattern,
                markersize=MARKERSIZE,
                color=line_color,
                marker="o",
            )
            ax[get_subplot_indices(plt_count, nrows, ncols)].set_title("# layers")
            if zcoord_type != "k":
                dxax = np.diff(
                    ax[get_subplot_indices(plt_count, nrows, ncols)].get_xlim()
                )[0]
                ax[get_subplot_indices(plt_count, nrows, ncols)].text(
                    1 + dxax * (0.05 + count * 0.26),
                    zcoord_m[-1],
                    "dz/2 = %d m" % (zm[-1] - ze[-1]),
                    horizontalalignment="left",
                    verticalalignment="center",
                    color=line_color,
                )
            if zlimit_type == "PBL":
                ax[get_subplot_indices(plt_count, nrows, ncols)].set_xlim(0, 35)

        if zcoord_type in ("p", "k"):
            ax[0, 0].invert_yaxis()
        if do_ylog:
            ax[0, 0].set_yscale("log")

        for r in range(ax.shape[0]):
            ax[r, 0].set_ylabel(zcoord_labels[zcoord_type])
        if not zlimit_type is None:
            ax[0, 0].set_ylim(zlimits[zcoord_type][zlimit_type])

        # plt.suptitle('%s (%d): p_s = %d hPa, z_s = %d m' % (levname, npz, ps/100., zs));
        # plt.suptitle('%s (%d)' % (levname, npz));
    fig.legend(
        *ax[0, 0].get_legend_handles_labels(),
        loc="upper center",
        bbox_to_anchor=(0.5, 0),
        ncol=len(ab_dict),
    )
    if save:
        plt_name = (
            "__".join([str.replace(" ", "_") for str in list(ab_dict.keys())])
            + "_levels.png"
        )
        plt.savefig(plt_name, dpi=DPI_FOR_PLOT)


def set_plotting_parameters():
    large = 24
    med = 20
    small = 16
    params = {
        "axes.titlesize": large,
        "legend.fontsize": med,
        "figure.figsize": (8, 4),
        "axes.labelsize": med,
        "axes.titlesize": large,
        "xtick.labelsize": small,
        "ytick.labelsize": small,
        "figure.titlesize": large,
        "axes.titlepad": 6,
    }
    plt.rcParams.update(params)
