# trace generated using paraview version 5.11.1
# import paraview
# paraview.compatibility.major = 5
# paraview.compatibility.minor = 11
from math import *
import h5py

# import the simple module from the paraview
from paraview.simple import *
from pathlib import Path

# disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()


work_dir = Path.cwd()
# work_dir = Path.home() / 'Desktop/Ke100_Pin3.0um'
# work_dir = Path('/Users/mnt/ceph/users/alamson/DATA/Chromatin/DynCondPaper/24-01-22_aLc1_line1600_Pin5um_Ke30_ks100_2patch')

assert work_dir.exists()
# assert analysis_dir.exists()

# with h5py.File(analysis_dir / 'raw_data.h5', 'r')  as h5d, (analysis_dir / 'test.txt').open('w') as file1:
#     # Examine the RunConfig file
#     # print( "###### RunConfig file string ########")
#     run_config_str = h5d.attrs['RunConfig']
#     # print(run_config_str, "\n")
#     file1.write(run_config_str)


# create a new 'PVD Reader'
sylinderpvtppvd = PVDReader(registrationName='Sylinderpvtp.pvd',
                            FileName=str(work_dir / 'result/Sylinderpvtp.pvd'))
sylinderpvtppvd.CellArrays = ['gid', 'group', 'isImmovable', 'radius', 'radiusCollision', 'length', 'lengthCollision', 'vel', 'omega', 'velCollision', 'omegaCollision', 'velBilateral', 'omegaBilateral',
                              'velNonBrown', 'omegaNonBrown', 'force', 'torque', 'forceCollision', 'torqueCollision', 'forceBilateral', 'torqueBilateral', 'forceNonBrown', 'torqueNonBrown', 'velBrown', 'omegaBrown', 'xnorm', 'znorm']
sylinderpvtppvd.PointArrays = ['endLabel']

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# get the material library
materialLibrary1 = GetMaterialLibrary()

# get display properties
sylinderpvtppvdDisplay = GetDisplayProperties(
    sylinderpvtppvd, view=renderView1)

# renderView1.ResetCamera(True)

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# create a new 'Legacy VTK Reader'
simBoxvtk = LegacyVTKReader(registrationName='simBox.vtk', FileNames=[
                            str(work_dir / 'result/simBox.vtk')])

# create a new 'PVD Reader'
proteinpvtppvd = PVDReader(registrationName='Proteinpvtp.pvd',
                           FileName=str(work_dir / 'result/Proteinpvtp.pvd'))
proteinpvtppvd.CellArrays = ['gid', 'tag']
proteinpvtppvd.PointArrays = ['idBind']

# create a new 'PVD Reader'
conBlockpvtppvd = PVDReader(registrationName='ConBlockpvtp.pvd',
                            FileName=str(work_dir / 'result/ConBlockpvtp.pvd'))
conBlockpvtppvd.CellArrays = [
    'oneSide', 'bilateral', 'delta0', 'gamma', 'kappa', 'Stress']
conBlockpvtppvd.PointArrays = ['gid', 'globalIndex', 'posIJ', 'normIJ']

# show data in view
proteinpvtppvdDisplay = Show(
    proteinpvtppvd, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
proteinpvtppvdDisplay.Representation = 'Surface'
proteinpvtppvdDisplay.ColorArrayName = [None, '']
proteinpvtppvdDisplay.SelectTCoordArray = 'None'
proteinpvtppvdDisplay.SelectNormalArray = 'None'
proteinpvtppvdDisplay.SelectTangentArray = 'None'
proteinpvtppvdDisplay.OSPRayScaleArray = 'idBind'
proteinpvtppvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
proteinpvtppvdDisplay.SelectOrientationVectors = 'None'
proteinpvtppvdDisplay.ScaleFactor = 0.049728404636776295
proteinpvtppvdDisplay.SelectScaleArray = 'None'
proteinpvtppvdDisplay.GlyphType = 'Arrow'
proteinpvtppvdDisplay.GlyphTableIndexArray = 'None'
proteinpvtppvdDisplay.GaussianRadius = 0.0024864202318388147
proteinpvtppvdDisplay.SetScaleArray = ['POINTS', 'idBind']
proteinpvtppvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
proteinpvtppvdDisplay.OpacityArray = ['POINTS', 'idBind']
proteinpvtppvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
proteinpvtppvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
proteinpvtppvdDisplay.PolarAxes = 'PolarAxesRepresentation'
proteinpvtppvdDisplay.SelectInputVectors = [None, '']
proteinpvtppvdDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
proteinpvtppvdDisplay.OSPRayScaleFunction.Points = [
    0.0, 0.0, 0.5, 0.0, 920.3745727539062, 0.30000001192092896, 0.5, 0.0, 4095.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
proteinpvtppvdDisplay.ScaleTransferFunction.Points = [
    -1.0, 0.0, 0.5, 0.0, 43.951139084439866, 0.30000001192092896, 0.5, 0.0, 199.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
proteinpvtppvdDisplay.OpacityTransferFunction.Points = [
    -1.0, 0.0, 0.5, 0.0, 43.951139084439866, 0.30000001192092896, 0.5, 0.0, 199.0, 1.0, 0.5, 0.0]

# show data in view
conBlockpvtppvdDisplay = Show(
    conBlockpvtppvd, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
conBlockpvtppvdDisplay.Representation = 'Surface'
conBlockpvtppvdDisplay.ColorArrayName = [None, '']
conBlockpvtppvdDisplay.SelectTCoordArray = 'None'
conBlockpvtppvdDisplay.SelectNormalArray = 'None'
conBlockpvtppvdDisplay.SelectTangentArray = 'None'
conBlockpvtppvdDisplay.OSPRayScaleArray = 'gid'
conBlockpvtppvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
conBlockpvtppvdDisplay.SelectOrientationVectors = 'None'
conBlockpvtppvdDisplay.ScaleFactor = 0.049728404636776295
conBlockpvtppvdDisplay.SelectScaleArray = 'None'
conBlockpvtppvdDisplay.GlyphType = 'Arrow'
conBlockpvtppvdDisplay.GlyphTableIndexArray = 'None'
conBlockpvtppvdDisplay.GaussianRadius = 0.0024864202318388147
conBlockpvtppvdDisplay.SetScaleArray = ['POINTS', 'gid']
conBlockpvtppvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
conBlockpvtppvdDisplay.OpacityArray = ['POINTS', 'gid']
conBlockpvtppvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
conBlockpvtppvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
conBlockpvtppvdDisplay.PolarAxes = 'PolarAxesRepresentation'
conBlockpvtppvdDisplay.SelectInputVectors = ['POINTS', 'normIJ']
conBlockpvtppvdDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
conBlockpvtppvdDisplay.OSPRayScaleFunction.Points = [
    0.0, 0.0, 0.5, 0.0, 920.3745727539062, 0.30000001192092896, 0.5, 0.0, 4095.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
conBlockpvtppvdDisplay.ScaleTransferFunction.Points = [
    0.0, 0.0, 0.5, 0.0, 44.72638338901767, 0.30000001192092896, 0.5, 0.0, 199.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
conBlockpvtppvdDisplay.OpacityTransferFunction.Points = [
    0.0, 0.0, 0.5, 0.0, 44.72638338901767, 0.30000001192092896, 0.5, 0.0, 199.0, 1.0, 0.5, 0.0]

# show data in view
simBoxvtkDisplay = Show(simBoxvtk, renderView1, 'UniformGridRepresentation')

# trace defaults for the display properties.
simBoxvtkDisplay.Representation = 'Outline'
simBoxvtkDisplay.ColorArrayName = ['POINTS', '']
simBoxvtkDisplay.SelectTCoordArray = 'None'
simBoxvtkDisplay.SelectNormalArray = 'None'
simBoxvtkDisplay.SelectTangentArray = 'None'
simBoxvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
simBoxvtkDisplay.SelectOrientationVectors = 'None'
simBoxvtkDisplay.SelectScaleArray = 'None'
simBoxvtkDisplay.GlyphType = 'Arrow'
simBoxvtkDisplay.GlyphTableIndexArray = 'None'
simBoxvtkDisplay.GaussianRadius = 0.05
simBoxvtkDisplay.SetScaleArray = [None, '']
simBoxvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
simBoxvtkDisplay.OpacityArray = [None, '']
simBoxvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
simBoxvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
simBoxvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
simBoxvtkDisplay.ScalarOpacityUnitDistance = 17.320508075688775
simBoxvtkDisplay.OpacityArrayName = [None, '']
simBoxvtkDisplay.ColorArray2Name = [None, '']
simBoxvtkDisplay.SliceFunction = 'Plane'
simBoxvtkDisplay.SelectInputVectors = [None, '']
simBoxvtkDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
simBoxvtkDisplay.OSPRayScaleFunction.Points = [
    0.0, 0.0, 0.5, 0.0, 920.3745727539062, 0.30000001192092896, 0.5, 0.0, 4095.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
simBoxvtkDisplay.ScaleTransferFunction.Points = [
    0.0, 0.0, 0.5, 0.0, 0.22475569542219934, 0.30000001192092896, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
simBoxvtkDisplay.OpacityTransferFunction.Points = [
    0.0, 0.0, 0.5, 0.0, 0.22475569542219934, 0.30000001192092896, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(sylinderpvtppvd)

# create a new 'Glyph'
glyph1 = Glyph(registrationName='Glyph1', Input=sylinderpvtppvd,
               GlyphType='Arrow')
glyph1.OrientationArray = ['POINTS', 'No orientation array']
glyph1.ScaleArray = ['POINTS', 'No scale array']
glyph1.ScaleFactor = 0.04972840463701927
glyph1.GlyphTransform = 'Transform2'

# show data in view
glyph1Display = Show(glyph1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
glyph1Display.ColorArrayName = [None, '']
glyph1Display.SelectTCoordArray = 'None'
glyph1Display.SelectNormalArray = 'None'
glyph1Display.SelectTangentArray = 'None'
glyph1Display.OSPRayScaleArray = 'endLabel'
glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph1Display.SelectOrientationVectors = 'None'
glyph1Display.ScaleFactor = 0.054452234692871575
glyph1Display.SelectScaleArray = 'None'
glyph1Display.GlyphType = 'Arrow'
glyph1Display.GlyphTableIndexArray = 'None'
glyph1Display.GaussianRadius = 0.0027226117346435784
glyph1Display.SetScaleArray = ['POINTS', 'endLabel']
glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph1Display.OpacityArray = ['POINTS', 'endLabel']
glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
glyph1Display.PolarAxes = 'PolarAxesRepresentation'
glyph1Display.SelectInputVectors = [None, '']
glyph1Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
glyph1Display.OSPRayScaleFunction.Points = [
    0.0, 0.0, 0.5, 0.0, 920.3745727539062, 0.30000001192092896, 0.5, 0.0, 4095.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
glyph1Display.ScaleTransferFunction.Points = [
    0.0, 0.0, 0.5, 0.0, 0.22475569542219934, 0.30000001192092896, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
glyph1Display.OpacityTransferFunction.Points = [
    0.0, 0.0, 0.5, 0.0, 0.22475569542219934, 0.30000001192092896, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on glyph1
glyph1.GlyphType = 'Sphere'
glyph1.GlyphType.ThetaResolution = 16
glyph1.GlyphType.PhiResolution = 16

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on glyph1
glyph1.GlyphMode = 'All Points'

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on glyph1
glyph1.ScaleArray = ['CELLS', 'radius']

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on glyph1
glyph1.ScaleFactor = 1.0

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on glyph1
glyph1.ScaleFactor = 2.0

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(glyph1Display, ('POINTS', 'gid'))

# rescale color and/or opacity maps used to include current data range
glyph1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# get 2D transfer function for 'gid'
gidTF2D = GetTransferFunction2D('gid')

# get color transfer function/color map for 'gid'
gidLUT = GetColorTransferFunction('gid')
gidLUT.TransferFunction2D = gidTF2D
gidLUT.RGBPoints = [0.0, 0.0, 0.0, 0.8, 45.04835164835165, 0.0, 0.0,
                    0.8, 153.95164835164834, 0.8, 0.3, 0.3, 1599.0, 0.0, 0.695201, 0.5]
gidLUT.ColorSpace = 'Step'
gidLUT.NanColor = [0.803922, 0.0, 0.803922]
gidLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'gid'
gidPWF = GetOpacityTransferFunction('gid')
gidPWF.Points = [0.0, 0.0, 0.5, 0.0, 44.72638338901767,
                 0.30000001192092896, 0.5, 0.0, 1599.0, 1.0, 0.5, 0.0]
gidPWF.ScalarRangeInitialized = 1

# Apply a preset using its name. Note this may not work as expected when
# presets have duplicate names.
gidLUT.ApplyPreset('Rainbow Uniform', True)
gidLUT.InvertTransferFunction()

# hide color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, False)

# rename source object
RenameSource('GID', glyph1)

gID = FindSource('GID')

# create a new 'Annotate Time Filter'
annotateTimeFilter1 = AnnotateTimeFilter(
    registrationName='AnnotateTimeFilter1', Input=gID)
annotateTimeFilter1Display = Show(
    annotateTimeFilter1, renderView1, 'TextSourceRepresentation')
annotateTimeFilter1.Scale = 0.5  # TODO make this something found from the data
# TODO make this hr:min:sec filter
annotateTimeFilter1.Format = 'Time: {time:.1f}'
annotateTimeFilter1Display.FontSize = 64
annotateTimeFilter1Display.WindowLocation = 'Upper Center'


# hide data in view
Hide(proteinpvtppvd, renderView1)

# hide data in view
Hide(conBlockpvtppvd, renderView1)

# hide data in view
Hide(simBoxvtk, renderView1)

# hide data in view
Hide(sylinderpvtppvd, renderView1)

# ================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
# ================================================================

# get layout
layout1 = GetLayout()

# split cell
layout1.SplitVertical(0, 0.5)

# set active view
SetActiveView(None)
# Adjust camera

# current camera placement for renderView1
# renderView1.CameraPosition = [0.24874988563699163,
# 0.00019760412025856405, 3.1148423714909423]
# renderView1.CameraFocalPoint = [
# 0.24874988563699163, 0.00019760412025856405, 4.898482660470327e-08]
# renderView1.CameraParallelScale = 0.28854242535090935

# Create a new 'Python View'
pythonView1 = CreateView('PythonView')
pythonView1.Script = """from paraview import python_view

def setup_data(view):
    pass

def render(view, width, height):
    figure = python_view.matplotlib_figure(width, height)
    ax = figure.add_subplot(1,1,1)
"""

# assign view to a particular cell in the layout
AssignViewToLayout(view=pythonView1, layout=layout1, hint=2)


# set active source
SetActiveSource(sylinderpvtppvd)

# show data in view
gIDDisplay = Show(sylinderpvtppvd, pythonView1, 'PythonRepresentation')

# find source
proteinpvtppvd = FindSource('Proteinpvtp.pvd')

# find source
conBlockpvtppvd = FindSource('ConBlockpvtp.pvd')

# find source
simBoxvtk = FindSource('simBox.vtk')

# find source
sylinderpvtppvd = FindSource('Sylinderpvtp.pvd')

# update the view to ensure updated data information
pythonView1.Update()

# Properties modified on pythonView1
pythonView1.Script = """from paraview import python_view
from paraview.simple import GetActiveView, GetAnimationScene, ExtractTimeSteps, GetTimeKeeper
from pathlib import Path
from paraview.numpy_support import vtk_to_numpy
import matplotlib.pyplot as plt
import h5py
import numpy as np
from matplotlib.collections import LineCollection


def gauss_weighted_contact(sep_mat, sigma=.010):
    return -np.power(sep_mat, 2) / (2. * (sigma * sigma))


def draw_vert_rainbow_line(ax, t, n_beads, cmap=\'jet_r\', lw=10):
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
                      contact_type="", vmax=25, label_flag=True):
    y = np.arange(contact_kymo.shape[0] + 1)
    # Add extra time point
    x = np.append(time_arr, [time_arr[-1] + time_arr[2] - time_arr[1]])
    X, Y = np.meshgrid(x, y)
    if contact_type == "log":
        c = ax.pcolorfast(X, Y, np.log(contact_kymo))
        if label_flag:
            _ = fig.colorbar(c, ax=ax, label="Log sum contact \
 probability")
    else:
        c = ax.pcolorfast(X, Y, contact_kymo, vmax=vmax)
        if label_flag:
            _ = fig.colorbar(c, ax=ax, label=r"Contact probability")
    if label_flag:
        _ = ax.set_title("Contact probabilty \'kymograph\'")
        ax.set_xlabel("Time $t$ (sec)")
        ax.set_ylabel("Bead index")
    
    ax.invert_yaxis()
    return


def setup_data(view):
    alens_stl = {{
        "axes.titlesize": 48,
        "axes.labelsize": 48,
        "lines.linewidth": 3,
        "lines.markersize": 10,
        "xtick.labelsize": 48,
        "ytick.labelsize": 48,
        "font.size": 48,
        "font.sans-serif": \'helvetica\',
        "text.usetex": False,
        \'mathtext.fontset\': \'cm\',
    }}
    plt.style.use(alens_stl)
    plt.rcParams[\'image.cmap\'] = \'YlOrRd\'

    pass


def render(view, width, height):
    fig = python_view.matplotlib_figure(width, height)
    axarr = []
    axarr += [fig.add_subplot(1, 2, 1)]
    axarr += [fig.add_subplot(1, 2, 2)]

    # fig, axarr = plt.subplots(1, 2, figsize=(20, 8))
    axarr[0].set_aspect(\'equal\')
    axarr[0].set_xlabel(r"Bead index")
    axarr[0].set_ylabel(r"Bead index")

    data_object = view.GetVisibleDataObjectForRendering(0)
    x = data_object.GetPoints().GetData()  # Get positions of points
    pos_arr = vtk_to_numpy(x)
    com_arr = .5*(pos_arr[:-1:2] + pos_arr[1::2])

    sep_mat = np.linalg.norm(
        com_arr[:, np.newaxis, :] - com_arr[np.newaxis, :, :], axis=2)
    contact_map = gauss_weighted_contact(sep_mat, sigma=.010)

    c = axarr[0].pcolorfast(contact_map, cmap=\'YlOrRd\', vmin=-25, vmax=0)
    axarr[0].invert_yaxis()
    fig.colorbar(c, ax=axarr[0],
                 # label=r"$\\log$(Inferred contact map) $\\sim$
                 # ($r_{{ij}}^2/2\\sigma^2$)")
                 label=r"Log contact probability $\\sim$ ($-r_{{ij}}^2$)",
                 )

    contact_path = Path("{}")
    assert contact_path.exists(), "Contact analysis file does not exist"
    with h5py.File(contact_path, \'r\') as h5d:
        CONTACT_KYMO = h5d[\'contact_kymo\'][...]
        TIME_ARR = h5d[\'time\'][...]

    # Get the current frame index
    current_frame_index = int(GetActiveView().ViewTime)
    print("Current timestep: ", current_frame_index)
    # c = axarr[1].pcolor(sep_mat, cmap=\'YlOrRd\')
    # axarr[1].set_title("Time {{:.2f}} sec".format(current_frame_index))
    plot_contact_kymo(fig, axarr[1], TIME_ARR, CONTACT_KYMO)
    # axarr[1].axvline(TIME_ARR[current_frame_index],
    #  color=\'w\', linestyle=\'--\')

    draw_vert_rainbow_line(
        axarr[1], TIME_ARR[current_frame_index], CONTACT_KYMO.shape[0])
    return python_view.figure_to_image(fig)
""".format(work_dir / 'analysis/contact_analysis.h5')

# ================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
# ================================================================

# --------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(4486, 2187)

# -----------------------------------
# saving camera placements for views

# current camera placement for renderView1
# renderView1.CameraPosition = [0.0, 0.0, 2.]
# renderView1.CameraFocalPoint = [
#     0.0, 0.0, 4.898482660470327e-08]
# renderView1.CameraParallelScale = 0.28854242535090935

# run the pipeline here to get the bounds
# rep = Show(bovr)
# rep.Representation = 'Outline'
# Render()

# --------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
# get layout
layout1 = GetLayout()

animationScene1 = GetAnimationScene()
animationScene1.AnimationTime = 0.0

# Enter preview mode
layout1.PreviewMode = [3840, 2160]
layout1.SetSize(3840, 2159)

renderView1.ResetCamera(True)

SaveAnimation(
#     str(work_dir / 'analysis/pv_flexible_filament_state.avi'), layout1, FrameRate=3, FrameWindow=[1, 3])
str(work_dir / 'analysis/pv_flexible_filament_state.avi'), layout1,
FrameRate=30, FrameWindow=[1, 1200])