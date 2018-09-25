import numpy as np
from uniform_bspline import UniformBspline

global pts
pts = [ [1,2,3,4,5,6,7,8,9],
        [2,3,5,7,3,6,8,9,1]]

def drawBSpline():
    '''
    Return the whole bspline and control point as bx, by, cx, cy
    Also return the velocity and acc and time as vx, vy, ax, ay, time
    '''
    global pts
    bspline = UniformBspline(pts)
    rg = bspline.getRegion()
    u = rg[0]

    while u<rg[1]:
        val = bspline.evaluate(u,0)

        val = bspline.evaluate(u,1)

        val = bspline.evaluate(u,2)

        u += 0.2

        print '----------------------------------------'

drawBSpline()