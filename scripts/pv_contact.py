from paraview import python_view
from paraview.simple import GetActiveView, GetAnimationScene, ExtractTimeSteps
from pathlib import Path
from paraview.numpy_support import vtk_to_numpy
import matplotlib.pyplot as plt
import h5py
import numpy as np
from matplotlib.collections import LineCollection


def gauss_weighted_contact(sep_mat, sigma=.010):
    return np.exp(-np.power(sep_mat, 2) / (2. * (sigma * sigma)))


def draw_vert_rainbow_line(ax, t, n_beads, cmap='jet', lw=10):
    c_arr = np.linspace(0, 1, n_beads)
    t_arr = np.ones(n_beads)*t

    points = np.array([t_arr, c_arr*n_beads]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap=cmap, norm=plt.Normalize(0, 1))
    lc.set_array(c_arr)
    lc.set_linewidth(lw)
    line = ax.add_collection(lc)
    lc.set_zorder(10)


def plot_contact_kymo(fig, ax, time_arr, contact_kymo,
                      contact_type="", vmax=10, label_flag=True):
    y = np.arange(contact_kymo.shape[0] + 1)
    # Add extra time point
    x = np.append(time_arr, [time_arr[-1] + time_arr[2] - time_arr[1]])
    X, Y = np.meshgrid(x, y)
    if contact_type == "log":
        c = ax.pcolorfast(X, Y, np.log(contact_kymo))
        if label_flag:
            _ = fig.colorbar(c, ax=ax, label="Log sum contact \n probability")
    else:
        c = ax.pcolorfast(X, Y, contact_kymo, vmax=vmax)
        if label_flag:
            _ = fig.colorbar(c, ax=ax, label=r"Contact probability")
    if label_flag:
        _ = ax.set_title("Contact probabilty 'kymograph'")
        ax.set_xlabel("Time $t$ (sec)")
        ax.set_ylabel("Bead index")
    return


def setup_data(view):
    alens_stl = {
        "axes.titlesize": 30,
        "axes.labelsize": 32,
        "lines.linewidth": 3,
        "lines.markersize": 10,
        "xtick.labelsize": 32,
        "ytick.labelsize": 32,
        "font.size": 32,
        "font.sans-serif": 'helvetica',
        "text.usetex": False,
        'mathtext.fontset': 'cm',
    }
    plt.style.use(alens_stl)
    plt.rcParams['image.cmap'] = 'YlOrRd'

    pass


def render(view, width, height):
    fig = python_view.matplotlib_figure(width, height)
    axarr = []
    axarr += [fig.add_subplot(1, 2, 1)]
    axarr += [fig.add_subplot(1, 2, 2)]

    # fig, axarr = plt.subplots(1, 2, figsize=(20, 8))
    axarr[0].set_aspect('equal')
    axarr[0].set_xlabel(r"Bead index")
    axarr[0].set_ylabel(r"Bead index")

    data_object = view.GetVisibleDataObjectForRendering(0)
    x = data_object.GetPoints().GetData()  # Get positions of points
    pos_arr = vtk_to_numpy(x)
    com_arr = .5*(pos_arr[:-1:2] + pos_arr[1::2])

    sep_mat = np.linalg.norm(
        com_arr[:, np.newaxis, :] - com_arr[np.newaxis, :, :], axis=2)
    contact_map = gauss_weighted_contact(sep_mat, sigma=.010)

    c = axarr[0].pcolorfast(contact_map, cmap='YlOrRd')
    fig.colorbar(c, ax=axarr[0],
                 # label=r"$\log$(Inferred contact map) $\sim$
                 # ($r_{ij}^2/2\sigma^2$)")
                 label=r"Log contact probability $\sim$ ($-r_{ij}^2$)")

    contact_path = Path.home() / \
        'projects/DATA/my_alens_data/NewStickyFlexibleFilament/analysis/contact_analysis.h5'
    assert contact_path.exists(), "Contact analysis file does not exist"
    with h5py.File(contact_path, 'r') as h5d:
        CONTACT_KYMO = h5d['contact_kymo'][...]
        TIME_ARR = h5d['time'][...]

    # Get the current frame index
    current_frame_index = int(GetActiveView().ViewTime)
    # c = axarr[1].pcolor(sep_mat, cmap='YlOrRd')
    # axarr[1].set_title("Time {:.2f} sec".format(current_frame_index))
    plot_contact_kymo(fig, axarr[1], TIME_ARR, CONTACT_KYMO)
    # axarr[1].axvline(TIME_ARR[current_frame_index],
    #  color='w', linestyle='--')

    draw_vert_rainbow_line(
        axarr[1], TIME_ARR[current_frame_index], CONTACT_KYMO.shape[0])
    return python_view.figure_to_image(fig)
