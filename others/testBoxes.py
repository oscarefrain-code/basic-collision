#============================================================================
#
#    testBoxes.py
#    Test the collision detector for two boxes. Visualize the contact points
#    in robotviewer to validate the collision detection.
#                                        Oscar E. Ramos Ponce, LAAS-CNRS
#
#============================================================================

import sys
sys.path.append('..')

from dynamic_graph.sot.core import *
from dynamic_graph.sot.dyninv import *
from dynamic_graph.script_shortcuts import optionalparentheses

from ThreadInterruptibleLoop import *
from time import sleep
import robotviewer


# -------------------------------------------------------------------
@loopInThread
def loop():
    inc()
runner=loop()

@optionalparentheses
def go(): runner.play()
@optionalparentheses
def stop(): runner.pause()
@optionalparentheses
def next(): inc()


# -------------------------------------------------------------------
# Initial position for the boxes
pbox1 = [0,0,0.1/2]; 
pbox2 = [0,0.05,0.1+0.024/2];

# Start the collision detector box to box
scene = BoxBoxCollisionDetector('scene')
scene.setTolerance(0.0001)
scene.addBox1((0.19, 0.115, 0.10),
              ((1.,0.,0.,pbox1[0]),(0.,1.,0.,pbox1[1]),(0.,0.,1.,pbox1[2]),(0.,0.,0.,1.)) )
scene.addBox2((0.23, 0.135, 0.024),
              ((1.,0.,0.,pbox2[0]),(0.,1.,0.,pbox2[1]),(0.,0.,1.,pbox2[2]),(0.,0.,0.,1.)) )

# Set the initial positions to robot-viewer
viewer = robotviewer.client('XML-RPC')
viewer.updateElementConfig('obstacle', [ pbox1[0], pbox1[1], pbox1[2], 0, 0, 0])
viewer.updateElementConfig('boxRF',    [ pbox2[0], pbox2[1], pbox2[2], 0, 0, 0])
# Hide all the contact points (under the ground)
MaxPt = 20
for i in range(MaxPt):
    name = 'point'+str(i+1)
    viewer.updateElementConfig(name, [0,0,-1,0,0,0])


# -------------------------------------------------------------------
def showVertices(Box, time):
    ''' To verify if the vertices are correctly calculated: VERIFIED for translation!!! '''
    scene.verticesB1.recompute(time)
    scene.verticesB2.recompute(time)
    if (Box==1):
        V = scene.verticesB1.value
    elif (Box==2):
        V = scene.verticesB2.value
    for i in range(8):
        name = 'point'+str(i+1)
        viewer.updateElementConfig(name, [V[i][0], V[i][1], V[i][2],0,0,0])    #Show the vertices
    for i in range(8,MaxPt):
        name = 'point'+str(i+1)
        viewer.updateElementConfig(name, [0,0,-1,0,0,0]) #Make dissappear the non-contact points


# Make one object move (translational motion)
def translateB2(dp, time):
    ''' dp = [delta_x, delta_y, delta_z] '''
    scene.transformationOut2.recompute(time)
    M = scene.transformationOut2.value 
    M = [list(M[0]), list(M[1]), list(M[2]), list(M[3])]
    M[0][3]+=dp[0]
    M[1][3]+=dp[1]
    M[2][3]+=dp[2]
    M = ( tuple(M[0]), tuple(M[1]), tuple(M[2]), tuple(M[3]) )
    scene.setTransformationB2(M)
    viewer.updateElementConfig('boxRF', [ M[0][3], M[1][3], M[2][3], 0, 0, 0])


def showCollision(time):
    scene.inContact.recompute(time)
    scene.contactPoints.recompute(time)
    if ( scene.inContact.value == 1): 
        Pc = scene.contactPoints.value                     # Contact points 'Pc' in the world frame 'w'
        for i in range(len(Pc)):
            name = 'point'+str(i+1)
            viewer.updateElementConfig(name, [Pc[i][0], Pc[i][1], Pc[i][2],0,0,0]) #Show the points
        for i in range(len(Pc),MaxPt):
            name = 'point'+str(i+1)
            viewer.updateElementConfig(name, [0,0,-1,0,0,0]) #Make dissappear the non-contact points

    
# -------------------------------------------------------------------
t = 0            
def inc():
    global t; t+=1
    sleep(0.005)
    #translateB2([0,0,-0.001], t)
    translateB2([0,0,-0.0005], t)
    showCollision(t)
    #showVertices(2, t)



