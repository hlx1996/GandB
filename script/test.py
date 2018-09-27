import numpy as np
from uniform_bspline import UniformBspline
from bspline import BSpline

global pts
pts = [ [1,2,3,4,5,6,7,8,9],
        [2,3,5,7,3,6,8,9,1]]
knots = [-4,-3,-2,-1,  0,2,4,6,8,10, 11,12,13,14]

def drawBSpline():
    '''
    Return the whole bspline and control point as bx, by, cx, cy
    Also return the velocity and acc and time as vx, vy, ax, ay, time
    '''
    global pts, knots
    bspline = BSpline(pts,knots)
    rg = bspline.getRegion()
#     print 'rg:', rg

    print bspline.evaluateByMat(1.0)
    print bspline.evaluateByMat(3.0)
    print bspline.evaluateByMat(5.0)
    print bspline.evaluateByMat(7.0)




drawBSpline()