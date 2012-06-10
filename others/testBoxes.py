#============================================================================
#
#    testBoxes.py
#    Test the collision detector for two boxes. Visualize the contact points
#    in robotviewer to validate the collision detection.
#                                        Oscar E. Ramos Ponce, LAAS-CNRS
#
#============================================================================

from dynamic_graph.sot.core import *
from dynamic_graph.sot.dyninv import *
#from dynamic_graph.script_shortcuts import optionalparentheses
#from dynamic_graph import plug
#from dynamic_graph.sot.dynamics import *


pbox1 = [0,0,0];     # Position for the first box
pbox2 = [0,0,0.4];   # Position for the second box 

# Start the collision detector box to box
scene = BoxBoxCollisionDetector('scene')
# scene.setTolerance(0.01)
# scene.addBox1((0.19, 0.115, 0.10),
#               ((1.,0.,0.,pbox1[0]),(0.,1.,0.,pbox1[1]),(0.,0.,1.,pbox1[2]),(0.,0.,0.,1.)) )
# scene.addBox2((0.23, 0.135, 0.024),
#               ((1.,0.,0.,pbox2[0]),(0.,1.,0.,pbox2[1]),(0.,0.,1.,pbox2[2]),(0.,0.,0.,1.)) )

# # Set the initial positions to robot-viewer
# robot.viewer.updateElementConfig('obstacle', [ pObst[0], pObst[1], pObst[2], 0, 0, 0])
# robot.viewer.updateElementConfig('boxRF',    [0,0,-1,0,0,0])      # View the right foot box

# # Change the position of the obstacle
# pObst = (0.03, -0.13, 0.03)
# robot.viewer.updateElementConfig('obstacle', [ pObst[0], pObst[1], pObst[2], 0, 0, 0])

# #rfCollision.setTolerance(0.05)
# wMar = matrix(dyn.rf.value)
# wMcomr = matrixToTuple(wMar*arMcomr)
# #rfCollision.addMovingBox((0.23, 0.135, 0.00001), (wMcomr))
# rfCollision.addBox2((0.23, 0.135, 0.024), (wMcomr))


# # Hide all the contact points (initialize them under the ground)
# for i in range(10):
#     name = 'point'+str(i+1)
#     robot.viewer.updateElementConfig(name, [0,0,-1,0,0,0]) #Show only the first point
# robot.viewer.updateElementConfig('boxRF',  [0,0,-1,0,0,0])      # View the right foot box

