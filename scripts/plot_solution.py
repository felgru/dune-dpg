#!/usr/bin/env pvpython

import sys
from paraview.simple import *

if len(sys.argv) != 3:
    print 'wrong number of command line arguments'
    print
    print 'Usage: plot_solution.py IN OUT'
    print
    print "Converts the 'theta' solution given by vtu file IN to a pdf file OUT"
    print

inputData = XMLUnstructuredGridReader(FileName=[sys.argv[1]])
inputData.PointArrayStatus = ['theta']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
renderView1.ViewSize = [1409, 807]
renderView1.UseOffscreenRendering = True
#renderView1.UseOffscreenRenderingForScreenshots = True

# display values of the solution on the third axis
solution = WarpByScalar(Input=inputData)
solution.Scalars = ['POINTS', 'theta']
solutionDisplay = Show(solution, renderView1)

# rescale color transfer function/color map
solutionLUT = GetColorTransferFunction('theta')
#solutionLUT.LoadPreset("Cool to Warm")
solutionLUT.RescaleTransferFunction(0.0, 1.0)
solutionLUT.LockScalarRange = 1
solutionPWF = GetOpacityTransferFunction('theta')
solutionPWF.RescaleTransferFunction(0.0, 1.0)

# show color bar/color legend
solutionDisplay.SetScalarBarVisibility(renderView1, True)

# set background color to white
renderView1.Background = [1.0, 1.0, 1.0]

# set color of text in the legend to black
solutionLUTColorBar = GetScalarBar(solutionLUT, renderView1)
solutionLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
solutionLUTColorBar.TitleFontFamily = 'Times'
solutionLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
solutionLUTColorBar.LabelFontFamily = 'Times'

# change title of color bar
solutionLUTColorBar.Title = ' '

# set camera properties
camera = GetActiveCamera()
camera.SetPosition((-1.17056, -2.26596, 1.73474))
camera.SetFocalPoint((0.5, 0.5, 0.5423))
camera.SetViewUp((0.214738, 0.274399, 0.937333))
camera.SetViewAngle(30)

Render()

# export view
ExportView(sys.argv[2], view=renderView1)
