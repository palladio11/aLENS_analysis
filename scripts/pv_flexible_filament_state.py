# trace generated using paraview version 5.11.1
# import paraview
# paraview.compatibility.major = 5
# paraview.compatibility.minor = 11

# import the simple module from the paraview
from paraview.simple import *
# disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'PVD Reader'
sylinderpvtppvd = PVDReader(registrationName='Sylinderpvtp.pvd',
                            FileName='/Users/alamson/projects/DATA/my_alens_data/NewStickyFlexibleFilament/result/Sylinderpvtp.pvd')
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

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# create a new 'Legacy VTK Reader'
simBoxvtk = LegacyVTKReader(registrationName='simBox.vtk', FileNames=[
                            '/Users/alamson/projects/DATA/my_alens_data/NewStickyFlexibleFilament/result/simBox.vtk'])

# create a new 'PVD Reader'
proteinpvtppvd = PVDReader(registrationName='Proteinpvtp.pvd',
                           FileName='/Users/alamson/projects/DATA/my_alens_data/NewStickyFlexibleFilament/result/Proteinpvtp.pvd')
proteinpvtppvd.CellArrays = ['gid', 'tag']
proteinpvtppvd.PointArrays = ['idBind']

# create a new 'PVD Reader'
conBlockpvtppvd = PVDReader(registrationName='ConBlockpvtp.pvd',
                            FileName='/Users/alamson/projects/DATA/my_alens_data/NewStickyFlexibleFilament/result/ConBlockpvtp.pvd')
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

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
renderView1.Update()

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
                    0.8, 153.95164835164834, 0.8, 0.3, 0.3, 199.0, 0.0, 0.695201, 0.5]
gidLUT.ColorSpace = 'Step'
gidLUT.NanColor = [0.803922, 0.0, 0.803922]
gidLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'gid'
gidPWF = GetOpacityTransferFunction('gid')
gidPWF.Points = [0.0, 0.0, 0.5, 0.0, 44.72638338901767,
                 0.30000001192092896, 0.5, 0.0, 199.0, 1.0, 0.5, 0.0]
gidPWF.ScalarRangeInitialized = 1

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
gidLUT.ApplyPreset('Rainbow Uniform', True)

# hide color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, False)

# rename source object
RenameSource('GID', glyph1)

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

# --------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(4486, 2414)

# -----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.CameraPosition = [0.24874988563699163,
                              0.00019760412025856405, 1.1148423714909423]
renderView1.CameraFocalPoint = [
    0.24874988563699163, 0.00019760412025856405, 4.898482660470327e-08]
renderView1.CameraParallelScale = 0.28854242535090935

# --------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
SaveScreenshot(
    '/Users/alamson/projects/DATA/my_alens_data/NewStickyFlexibleFilament/analysis/pv_flexible_filament_state.png')
