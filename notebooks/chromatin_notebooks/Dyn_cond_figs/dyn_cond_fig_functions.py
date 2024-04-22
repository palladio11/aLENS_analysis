# Basic useful imports
import re
import time
import yaml
from pprint import pprint
from pathlib import Path
import h5py
import warnings

# Data manipulation
import numpy as np
from scipy.special import erf
from scipy.integrate import quad
from scipy.optimize import root_scalar
import scipy.stats as stats
from scipy.signal import savgol_filter
from scipy.spatial import ConvexHull

# Visualization
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, RegularPolygon, FancyArrowPatch, ArrowStyle
from matplotlib.ticker import (
    MultipleLocator,
    FormatStrFormatter,
    AutoMinorLocator,
    NullFormatter,
)
import matplotlib.colors as mcolors

# Clustering stuff
from sklearn.cluster import MeanShift, estimate_bandwidth, DBSCAN, OPTICS
from itertools import cycle
# plt.cm.tab20.colors

# From alens_analysis.py
import alens_analysis as aa
import alens_analysis.chromatin as aac
import alens_analysis.chromatin.chrom_analysis as ca
import alens_analysis.chromatin.chrom_condensate_analysis as cca
import alens_analysis.chromatin.chrom_graph_funcs as cgf
from alens_analysis import cluster_analysis as cla

from alens_analysis.colormaps import register_cmaps

# Locations
ws_path = Path("/home/alamson/DATA/Chromatin/")
mnt_path = Path.home() / "projects/DATA/Chromatin/"
ceph_path = Path.home() / "ceph/DATA/Chromatin/"

graph_sty = {
    "axes.titlesize": 20,
    "axes.labelsize": 24,
    "lines.linewidth": 2,
    "lines.markersize": 2,
    "xtick.labelsize": 24,
    "ytick.labelsize": 24,
    "font.size": 20,
    "font.sans-serif": "Helvetica",
    "text.usetex": False,
    "mathtext.fontset": "cm",
}

cluster_similarity_threshold = 0.4
nskip = 10  # Time snapshot skips for cluster finding. = 5 secs
vmax = 40  # Max colorbar value in kymographs
tree_length = 30  # min length of a cluster tree in time snapshots. = 15 sec

colors = cycle(mcolors.XKCD_COLORS.keys())

register_cmaps()
# plt.rcParams['image.cmap'] = 'emct8'
# plt.rcParams['image.cmap'] = 'warm'
plt.rcParams["image.cmap"] = "YlOrRd"
# plt.rcParams['image.cmap'] = 'twilight'
# plt.rcParams['image.cmap'] = 'coolwarm'
# plt.rcParams['image.cmap'] = 'RdYlBu_r'


def plot_confidence_int(
    ax, time_arr, mean, std_dev, num_runs=12, color="b", ci=0.95, label="Mean"
):
    degrees_freedom = num_runs - 1
    confidence_interval = (
        stats.t.ppf((1 + ci) / 2.0, degrees_freedom) * std_dev / np.sqrt(num_runs)
    )

    _ = ax.plot(time_arr, mean, label=label, color=color)
    _ = ax.fill_between(
        time_arr,
        mean - confidence_interval,
        mean + confidence_interval,
        color=color,
        alpha=0.1,
    )


def cluster_analysis_graph(sim_path, part_min=40):
    flat_time_arr = []
    flat_clust_cent_arr = []
    flat_clust_ind_arr = []
    num_clusters_list = []
    num_cluster_beads_list = []

    # Length of chain and time to look at
    ss_ind = 1
    end_ind = None
    start_bead = 0
    end_bead = None
    with h5py.File(next(sim_path.glob("analysis/raw*.h5")), "r") as h5_data:
        time_arr = h5_data["time"][ss_ind:end_ind]
        print(time_arr.shape)
        sy_dat = h5_data["raw_data/sylinders"][start_bead:end_bead, :, ss_ind:end_ind]
        com_arr = 0.5 * (sy_dat[:, 2:5, :] + sy_dat[:, 5:8, :])

    h5_clust_file = sim_path / "analysis/cluster_analysis.h5"
    with h5py.File(h5_clust_file, "r") as h5_data:
        cluster_grp = h5_data["clusters"]
        time_arr = h5_data["time"][...]
        time_grp_list = sorted(cluster_grp.values(), key=lambda x: x.attrs["time"])
        clusters = []
        for tg in time_grp_list:
            clusters += [[cla.Cluster(h5_data=c) for c in tg.values()]]

    bead_ind_arr = np.zeros((com_arr.shape[0], com_arr.shape[2]))
    one_mask = np.ones(com_arr.shape[0])
    for c, clust_grp in enumerate(clusters):
        # Secondary thresholding
        clust_grp = [clust for clust in clust_grp if len(clust.part_ids) > part_min]
        num_clusters_list += [len(clust_grp)]
        num_beads = 0
        for i, clust in enumerate(clust_grp):
            flat_time_arr += [clust.time]
            flat_clust_cent_arr += [clust.center]
            flat_clust_ind_arr += [clust.part_ids]
            num_beads += len(clust.part_ids)
            bead_ind_arr[clust.part_ids, c] += one_mask[clust.part_ids]

        num_cluster_beads_list += [num_beads]

    # for c, clust_grp in enumerate(clusters):
    #     for clust in clust_grp:
    #         bead_ind_arr[clust.part_ids,c] += one_mask[clust.part_ids]
    # bead_ind_arr[:,:-1] += bead_ind_arr[:,1:]
    # bead_ind_arr *= .5

    fig, axarr = plt.subplots(2, 2, figsize=(20, 15))
    X, Y = np.meshgrid(time_arr, np.arange(com_arr.shape[0]))
    c = axarr[1, 0].pcolor(X, Y, bead_ind_arr, shading="nearest")

    flat_clust_cent_arr = np.asarray(flat_clust_cent_arr)
    _ = axarr[0, 0].plot(flat_time_arr, flat_clust_cent_arr[:, 0], ".")
    _ = axarr[0, 1].scatter(time_arr, num_clusters_list)
    _ = axarr[1, 1].scatter(time_arr, num_cluster_beads_list)
    _ = axarr[0, 0].set_ylabel("$x$-position ($\mu$m)")
    _ = axarr[0, 1].set_ylabel("Number of clusters")
    _ = axarr[1, 0].set_ylabel("Bead index")
    _ = axarr[1, 1].set_ylabel("Number of beads in clusters")
    for ax in axarr.flatten():
        _ = ax.set_xlabel("Time (sec)")
    fig.tight_layout()


def kymo_contact_graph(fig, ax, sim_path):
    h5_contact_file = sim_path / "analysis/contact_analysis.h5"

    with h5py.File(h5_contact_file, "r") as h5_data:
        time_arr = h5_data["time"][...]
        contact_kymo = h5_data["contact_kymo"][...]

    y = np.arange(contact_kymo.shape[0] + 1)
    x = np.append(time_arr, [time_arr[-1] + time_arr[2] - time_arr[1]])
    X, Y = np.meshgrid(x, y)
    c = ax.pcolorfast(X, Y, contact_kymo, vmax=40)
    _ = fig.colorbar(c, ax=ax, label=r"Contact number")
    # _ = x.set_xlabel("Time $t$ [sec]")
    # _ = x.set_ylabel("Bead index")
    _ = ax.set_ylim(0, 1600)

    _ = ax.invert_yaxis()


def kymo_and_cluster_graph(
    sim_path,
    cluster_similarity_threshold=0.4,
    nskip=10,
    vmax=40,
    tree_length=30,
):
    h5_contact_file = sim_path / "analysis/contact_analysis.h5"

    fig, axarr = plt.subplots(1, 1, figsize=(10, 6))
    _ = axarr.set_ylim(0, 1600)

    with h5py.File(h5_contact_file, "r") as h5_data:
        time_arr = h5_data["time"][...]
        contact_kymo = h5_data["contact_kymo"][...]

        cgf.plot_contact_kymo(fig, axarr, time_arr, contact_kymo, vmax=vmax)

    # Cluster analysis
    h5_clust_file = sim_path / "analysis/cluster_analysis.h5"
    with h5py.File(h5_clust_file, "r") as h5_data:
        cluster_grp = h5_data["clusters"]
        time_arr = h5_data["time"][...]
        time_grp_list = sorted(cluster_grp.values(), key=lambda x: x.attrs["time"])
        clusters = []
        for tg in time_grp_list:
            clusters += [[cla.Cluster(h5_data=c) for c in tg.values()]]
    root_clusters = cla.find_descendants(
        clusters, thresh=cluster_similarity_threshold, nskip=nskip
    )

    trees = []
    tree_id_gen = aa.helpers.gen_id()
    for root in root_clusters:
        tree = cla.ClusterTree(next(tree_id_gen))
        tree.add_recursive(root)
        trees += [tree]

    fig, axarr = plt.subplots(1, 2, figsize=(16, 6))

    # Graph all clusters
    for tree, color in zip(trees, colors):
        if len(tree.clusters) < tree_length:
            continue
        for clust in tree.clusters:
            _ = axarr[0].plot(
                [clust.time] * len(clust.part_ids),
                clust.part_ids,
                color=color,
                markersize=0.1,
                marker=".",
                linestyle="None",
            )

    # Graph largest branch
    biggest_tree = max(trees, key=lambda x: x.get_main_clust_branch()[0].mass_hist)
    biggest_tree.update_branch_roots()
    branch_roots = biggest_tree.get_branch_roots()

    for root, color in zip(branch_roots, colors):
        branch_clusters = root.get_largest_branch()
        for clust in branch_clusters:
            _ = axarr[1].plot(
                [clust.time] * len(clust.part_ids),
                clust.part_ids,
                color=color,
                markersize=0.1,
                marker=".",
                linestyle="None",
            )

    _ = axarr[0].set_ylabel("Bead index")
    for ax in axarr:
        _ = ax.set_ylim(0, 1600)
        _ = ax.set_xlim(0, time_arr[-1])
        _ = ax.set_xlabel("Time")


def calc_droplet_only_fe(ld, alpha, gamma, nu):
    return -(nu * alpha * ld) + (4.0 * np.pi * gamma) * np.power(
        (3.0 * alpha * ld) / (4.0 * np.pi), 2.0 / 3.0
    )


def calc_polymer_only_fe(ld, L, Lc, kappa):
    Lp = Lc - ld
    return 0.25 * kappa * ((Lp**2) / (Lp - L) + L**2 / Lp - (Lp + L))


def calc_polymer_only_fe_deriv(ld, L, Lc, kappa):
    Lp = Lc - ld
    return (
        (0.25 * kappa * L**2)
        * ((L**2) - (2.0 * Lp * L) + 2.0 * (Lp**2))
        / ((Lp**2) * ((L - Lp) ** 2))
    )


def calc_droplet_only_fe_deriv(ld, alpha, gamma, nu):
    return (
        2.0 * gamma * np.power(4.0 * np.pi * alpha * alpha / (3 * ld), 1.0 / 3.0)
    ) - (nu * alpha)


# Derivative of free energy of the system taken with respect to l1
def two_cond_free_energy_deriv(
    l1: float,
    l2: float,
    L: float,
    Lc: float,
    alpha: float,
    gamma: float,
    nu: float,
    kappa: float,
) -> float:
    """_summary_

    Parameters
    ----------
    l1 : float
        _description_
    l2 : float
        _description_
    L : float
        _description_
    Lc : float
        _description_
    alpha : float
        _description_
    gamma : float
        _description_
    nu : float
        _description_
    kappa : float
        _description_

    Returns
    -------
    float
        _description_
    """

    fe_deriv_droplet = calc_droplet_only_fe_deriv(l1, alpha, gamma, nu)
    fe_deriv_polymer = calc_polymer_only_fe_deriv(l1 + l2, L, Lc, kappa)
    return fe_deriv_droplet + fe_deriv_polymer


def free_energy_two_identical_conds_continuous_deriv(
    ld, L, Lc=1.0, nu=1.0, alpha=1.0, gamma=1.0, kappa=1.0
):
    """
    Free energy of system when you have two identical condensates with ld length of chain in each

    Parameters
    ----------
    ld : float
        amount of chain in a single blob
    L : float
        Separation of filament ends
    Lc : float, optional
        _description_, by default 1.0
    nu : float, optional
        _description_, by default 1.0
    alpha : float, optional
        _description_, by default 1.0
    gamma : float, optional
        _description_, by default 1.0
    kappa : float, optional
        _description_, by default 1.0

    Returns
    -------
    _type_
        _description_
    """

    return (1.0 / 6.0) * (
        # Bulk term
        -12.0 * alpha * nu
        # Polymer term
        + (3.0 * kappa * (L**2))
        * (L**2 - 2 * L * (Lc - 2 * ld) + 2 * ((Lc - 2 * ld) ** 2))
        / ((Lc - 2 * ld) ** 2 * (L - Lc + 2 * ld) ** 2)
        # Surface tension term
        + (8 * np.power(6, 2.0 / 3.0) * alpha * gamma * np.cbrt(np.pi / (alpha * ld)))
    )


def calc_max_length_in_two_condensates(
    L: float = 1.0,
    Lc: float = 1.0,
    nu: float = 1.0,
    alpha: float = 1.0,
    gamma: float = 1.0,
    kappa: float = 1.0,
    **kwargs,
):
    epsilon = 0.000001
    ld_lower_bound = epsilon  # Can't be zero,
    ld_upper_bound = (
        (Lc - (L + epsilon)) * 0.5
    )  # Half because there are two condensates contributing to total length in droplet
    while (
        free_energy_two_identical_conds_continuous_deriv(
            ld_lower_bound, L, Lc, nu, alpha, gamma, kappa
        )
        > 0
    ):
        ld_lower_bound += 0.05

    if (
        ld_lower_bound > ld_upper_bound
        or free_energy_two_identical_conds_continuous_deriv(
            ld_upper_bound, L, Lc, nu, alpha, gamma, kappa
        )
        < 0
    ):
        # No max was found, condensates would not exist
        print("No max found with end separation", L)
        return 0

    result = root_scalar(
        free_energy_two_identical_conds_continuous_deriv,
        method="brentq",
        args=(L, Lc, nu, alpha, gamma, kappa),
        bracket=[ld_lower_bound, ld_upper_bound],
    )
    return result.root


def two_cond_size_continuous_deriv(t, state, nu, gamma, alpha, kappa, L, Lc, b, beta):
    l1, l2 = state

    dAdl1 = two_cond_free_energy_deriv(l1, l2, L, Lc, alpha, gamma, nu, kappa)
    dAdl2 = two_cond_free_energy_deriv(l2, l1, L, Lc, alpha, gamma, nu, kappa)
    dl1 = 2.0 * b * (np.exp(-beta * b * dAdl1) - 1.0)
    dl2 = 2.0 * b * (np.exp(-beta * b * dAdl2) - 1.0)
    return [dl1, dl2]


def convert_e_to_x10n(text):
    # Find numbers in scientific notation
    matches = re.findall(r"\b\d+\.\d+e[+-]\d+\b", text)
    for match in matches:
        # Split the number into the coefficient and the exponent
        coefficient, exponent = match.split("e")
        # Convert the number to 'x10^n' format
        x10n_format = rf"${coefficient}\times 10^{int(exponent)}$"
        # Replace the original number with the 'x10^n' format
        text = text.replace(match, x10n_format)
    return text


def tmerge_exact(x_sep, y_com, L, Lc, nu, gamma, alpha, kappa, b, beta, kmodes=100):
    # Find max number of beads in two condensates
    max_ld = calc_max_length_in_two_condensates(L, Lc, nu, alpha, gamma, kappa)

    # Get max separation because condensates take up beads
    max_sep = Lc - 2 * max_ld
    assert x_sep < max_sep, "Separation too large for condensate size"
    # Translate to center of mass
    y_com += 0.5 * (max_sep - Lc)
    assert y_com > 0 and y_com < max_sep, "Center of mass is no longer on the chain"

    # Convert to bead number
    max_bead_in_cond = int(max_ld / b)

    n_beads = int(Lc / b)

    # Find free energy of the two condensates
    full_free_energy = calc_free_energy_two_blobs(
        max_bead_in_cond, max_bead_in_cond, n_beads, L, alpha, gamma, nu, kappa, bd=b
    )
    free_energy_minus_bead = calc_free_energy_two_blobs(
        max_bead_in_cond,
        max_bead_in_cond - 1,
        n_beads,
        L,
        alpha,
        gamma,
        nu,
        kappa,
        bd=b,
    )
    energy_diff = full_free_energy - free_energy_minus_bead

    # Get diffusion constant (This seems to be off by a factor of 4)
    D = 0.5 * b * b  # 2x a single condensate

    # Rescale
    x = x_sep / (2.0 * max_sep)
    y = y_com / (max_sep)

    # First term
    T_o = ((max_sep**2) * x / D) * (1.0 - x)

    # Second term
    T_1 = 0
    for i in range(kmodes):
        k = 2 * i + 1
        T_1 += (
            (np.sin(k * np.pi * x) / (k**3))
            * (np.sinh(k * y * np.pi) + np.sinh(k * (1 - y) * np.pi))
            / np.sinh(k * np.pi)
        )

    return T_o - ((8.0 * (max_sep**2) / (D * (np.pi**3))) * T_1)


FE_PREFACT = np.power(36.0 * np.pi, 1.0 / 3.0)


def calc_free_energy_two_blobs(
    l1_beads: int,
    l2_beads: int,
    n_beads: int,
    L: float = 1.0,
    alpha: float = 1.0,
    gamma: float = 1.0,
    nu: float = 1.0,
    kappa: float = 1.0,
    bd: float = 1.0,
    **kwargs,
) -> float:
    """Calculate and return the total free energy for the current configuration

    Parameters
    ----------
    n_beads : int
        Contour length of polymer in bead number
    l1_beads : int
        Length of chain in blob 1  in bead number
    l2_beads : int
        Length of chain in blob 2  in bead number
    L : float, optional
        Filament end separation [Length], by default 1
    alpha : float, optional
        Condensate packing volume per length of chain [Length^2], by default 1
    gamma : float, optional
        Surface tension [Force/Length], by default 1
    nu : float, optional
        Free energy per volume of condensate [Force/Length^2], by default 1
    kappa : float, optional
        Filament flexibility [Force], by default 1
    bd : float, optional
        Diameter of a bead [Length], by default 1

    Returns
    -------
    float
        _description_
    """

    ltot_beads = l1_beads + l2_beads
    Lp = (n_beads - ltot_beads) * bd

    term1 = (FE_PREFACT * gamma) * (
        np.power(l1_beads * bd * alpha, 2.0 / 3.0)
        + np.power(l2_beads * bd * alpha, 2.0 / 3.0)
    )
    term2 = 0.25 * kappa * L * L * ((L - (2.0 * Lp)) / (Lp * (L - Lp)))
    term3 = ltot_beads * alpha * nu * bd
    return term1 + term2 - term3
