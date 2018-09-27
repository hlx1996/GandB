import numpy as np
import math


'''
    This script implement a b-spline class, which is initialize by some control point which is initialized by some control points
    It is a general class for b-spline, so the knots vector is not necessarily uniform.
    We also use the general matrix representation to compute the b-spline, instead of using the precomputed coefficient
'''

class BSpline:
    def __init__(self, pts):
        # p is predefined
        self.p = 4
        # pts is in the form of [[x....],[y.....]], copy the first and last point p time
        ptx, pty = pts[0], pts[1]
        cptx, cpty = [], []
        for i in range(self.p):
            cptx.append(ptx[0])
            cpty.append(pty[0])
        for i in range(len(ptx)):
            cptx.append(ptx[i])
            cpty.append(pty[i])
        for i in range(self.p):
            cptx.append(ptx[len(ptx)-1])
            cpty.append(pty[len(pty)-1])
        # print cptx, cpty

        self.control_points = [cptx, cpty]
        # The control point numebr should add 2*p since begin and end point repeat p+1 time
        self.n = len(pts[0])+2* self.p-1
        # compute m
        self.m = self.n + self.p + 1
        # create knots vector for the control point
        # here we want up =0, up+1 is close to up, um-p-1 close to um-p 
        knots = []
        for i in range(self.m + 1):
            if i <= self.p:
                u = 0 - (self.p - i)
            elif i > self.p and i <= (self.m - self.p):
                u = knots[-1] + 2
            elif i > (self.m -self.p):
                u = knots[-1] + 1
            knots.append(u)
        self.u = knots

        # print self.n, self.m
        # print self.control_points, len(self.control_points[0]) 
        # print self.u, len(self.u)

        # The M5 matrix used for evaluation
        self.m5 = np.matrix([[1.0,11.0,11.0,1.0,0],[-4.0,-12.0,12.0,4.0,0],[6.0,-6.0,-6.0,6.0,0],[-4.0,12.0,-12.0,4.0,0],[1.0,-4.0,6.0,-4.0,1.0]])/(4.0*3.0*2)
        # print self.m5

    def getRegion(self):
        rg = [self.u[self.p], self.u[self.m-self.p]]
        return rg


    def evaluateByMat(self, u):
        '''
        Here we use the matrix representation to compute the bspline
        Since self.p = 4, we should use M5
        we first decide which [ui, ui+1] u lay in, then shift it to [0,1],
        and use the p = u'*M5*p
        '''
        if u < self.u[self.p] or u > self.u[self.m-self.p]:
            return

        # determine which [ui,ui+1] lay in
        idx = self.p
        while True:
            if self.u[idx+1] >= u:
                break
            idx += 1
        # print idx

        # Then shift u to [0,1]
        u = (u-self.u[idx])/(self.u[idx+1]-self.u[idx])
        # print u

        # use p = u'*M5*p
        up = np.matrix([1,u,u**2,u**3,u**4])
        # print up

        pi = self.control_points
        px = np.matrix([pi[0][idx-self.p],pi[0][idx-self.p+1],pi[0][idx-self.p+2],pi[0][idx-self.p+3],pi[0][idx]])
        py = np.matrix([pi[1][idx-self.p],pi[1][idx-self.p+1],pi[1][idx-self.p+2],pi[1][idx-self.p+3],pi[1][idx]])

        val = [0.0,0.0]
        val[0] = float(up*self.m5*np.transpose(px))
        val[1] = float(up*self.m5*np.transpose(py))
        # print val
        return val
            