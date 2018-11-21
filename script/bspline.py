import numpy as np
import math


'''
    This script implement a b-spline class, which is initialize by some control point which is initialized by some control points and order p
    It is a general class for b-spline
    We also use the general matrix representation to compute the b-spline, instead of using the precomputed coefficient
    The derivative of the b-spline can also be computed, since the derivative is also an b-spline
    auto_extend is used when the first and last control points need to be repeated p times. This is not need when create derivative 
'''

class BSpline:
    def __init__(self, pts, order, auto_extend=True):
        # p is predefined
        self.p = order
        if auto_extend:            
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
        else:
            # no need to repeat the first and last control points, just copy
            self.control_points = pts
            self.n = len(pts[0]) - 1
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
        print self.u, len(self.u)

        # The M matrix used for evaluation
        m5 = np.matrix([[1.0,11.0,11.0,1.0,0],[-4.0,-12.0,12.0,4.0,0],[6.0,-6.0,-6.0,6.0,0],[-4.0,12.0,-12.0,4.0,0],[1.0,-4.0,6.0,-4.0,1.0]])/(4.0*3.0*2)
        m4 = np.matrix([[1.0,4.0,1.0,0.0],[-3.0,0.0,3.0,0.0],[3.0,-6.0,3.0,0.0],[-1.0,3.0,-3.0,1.0]])/(3.0*2.0)
        m3 = np.matrix([[1.0,1.0,0.0],[-2.0,2.0,0.0],[1.0,-2.0,1.0]])/(2.0)
        self.mdict = {3:m3, 4:m4, 5:m5}

    def getRegion(self):
        rg = [self.u[self.p], self.u[self.m-self.p]]
        return rg

    def getU(self,u):
        '''
        Input: u in [up, um-p]
        Output: [1,u,u^2,...,u^p]
        '''
        uv = []
        for i in range(self.p+1):
            element = u**i
            uv.append(element)

        return np.matrix(uv)

    def getPi(self,idx):
        '''
        Input: index of idx -> ui
        Output: [pi-p,pi-p+1,...,pi]
        '''
        pi = self.control_points
        px, py = [], []
        for i in range(self.p+1):
            elx = pi[0][idx-self.p+i]
            ely = pi[1][idx-self.p+i]
            px.append(elx)
            py.append(ely)

        return np.transpose(np.matrix(px)), np.transpose(np.matrix(py))


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
        up = self.getU(u)
        px, py = self.getPi(idx)
        # print up

        val = [0.0,0.0]
        val[0] = float(up*self.mdict[self.p+1]*px)
        val[1] = float(up*self.mdict[self.p+1]*py)
        # print val
        return val

    def getDerivative(self):
        '''
        The derivative of a b-spline is also a b-spline, its order become p-1
        control point Qi = p*(Pi+1-Pi)/(ui+p+1-ui+1)
        '''
        
        # first calculate the new control points
        ptx = self.control_points[0]
        pty = self.control_points[1]
        d1x, d1y = [], []
        for i in range(len(ptx)-1):
            x = self.p * (ptx[i+1]-ptx[i])/(self.u[i+self.p+1]-self.u[i+1])
            y = self.p * (pty[i+1]-pty[i])/(self.u[i+self.p+1]-self.u[i+1])
            d1x.append(x)
            d1y.append(y)

        d1 = [d1x, d1y]
        
        derivative = BSpline(d1, self.p-1, False)

        return derivative


            