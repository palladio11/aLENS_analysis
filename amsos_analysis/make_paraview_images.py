#!/usr/bin/env python

"""@package docstring
File: make_paraview_images.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""

import sys
from paraview.simple import (paraview, LoadState, GetAnimationScene,
                             GetTimeKeeper, FindViewOrCreate, GetLayoutByName,
                             SetActiveView, SaveAnimation)


def main(state_file):
    """TODO: Docstring for main.

    @param state_file TODO
    @return: TODO

    """
    # disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()
    # load state
    LoadState(state_file,
              LoadStateDataFileOptions='Search files under specified directory',
              DataDirectory='./result',
              OnlyUseFilesInDataDirectory=1,
              simBoxvtkFileNames=['./result/simBox.vtk'],
              SylinderpvtppvdFileName='./result/Sylinderpvtp.pvd',
              ProteinpvtppvdFileName='./result/Proteinpvtp.pvd',
              ConBlockpvtppvdFileName='./result/ConBlockpvtp.pvd')
    # get animation scene
    animationScene1 = GetAnimationScene()
    # get the time-keeper
    timeKeeper1 = GetTimeKeeper()
    # find view
    renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')
    # uncomment following to set a specific view size
    # renderView1.ViewSize = [3199, 810]
    # get layout
    layout1_1 = GetLayoutByName("Layout #1")
    # set active view
    SetActiveView(renderView1)
    # Enter preview mode
    layout1_1.PreviewMode = [1920, 1280]
    # find view
    # renderView2 = FindViewOrCreate('RenderView2', viewtype='RenderView')
    # uncomment the following to render all views
    #times = reader.TimestepValues
    #view.ViewTime = times[n]
    # RenderAllViews()
    # alternatively, if you want to write images, you can use
    # SaveScreenshot(...).
    SaveAnimation('./result/PNG/image.png',
                  SaveAllViews=1,
                  viewOrLayout=layout1_1,
                  SeparatorWidth=5,
                  #    FrameWindow=[100, 3300],
                  SuffixFormat=".%06d")


#########################################
if __name__ == "__main__":
    main(sys.argv[1])
