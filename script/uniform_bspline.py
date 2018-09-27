import numpy as np
import math
from coefficent import coe


'''
    This script implement a b-spline class, which is initialize by some control point
Its p is fixed, n and m will be computed accordingly.

    It should be able to return a list [u_p, u_m-p] indicating the effective region.

    Then given a 'u' between [u_p, u_m-p], it should return the value of that point, including position, velocity and acceleration
'''

class UniformBspline:
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
        print cptx, cpty

        self.control_points = [cptx, cpty]
        print len(self.control_points[0])
        # The control point numebr should add 2*p since begin and end point repeat p+1 time
        self.n = len(pts[0])+2* self.p-1
        # compute m
        self.m = self.n + self.p + 1
        # pre-computed coefficient of uniform b-spline
        self.coe = coe

    def getRegion(self):
        rg = [self.p, self.m-self.p]
        return rg

    def getTimeVector(self,u):
        uv = []
        for i in range(self.p+1):
            uv.append(u**(self.p-i))
        # print uv
        return uv

    def evaluate(self, u, order):
        '''
        decide whether u lie between [up, u_m-p]
        return the value in the form of [x, y]
        order 0 for pos, 1 for vel, 2 for acc
        '''
        if u < self.p or u > (self.m-self.p):
            return

        # compute value of u
        xu = 0.0
        yu = 0.0
        uv = self.getTimeVector(u)
        s = int(math.floor(u-self.p))

        for i in range(self.p+1):
            '''
            transform the coefficient if order is larger than 0, T matrix transform the
            coefficent of pos to vel
            '''
            tmat = np.matrix([[0,0,0,0,0],[4,0,0,0,0],[0,3,0,0,0],[0,0,2,0,0],[0,0,0,1,0]])
            nip = 0.0
            if order == 0:
                nip = uv * np.transpose(np.matrix(self.coe[s][i]))
            elif order == 1:
                nip = uv *tmat* np.transpose(np.matrix(self.coe[s][i]))
            elif order == 2:
                nip = uv *tmat*tmat* np.transpose(np.matrix(self.coe[s][i]))
            xu += self.control_points[0][s+i] * nip
            yu += self.control_points[1][s+i] * nip
        val = [float(xu),float(yu)]

        return val


def main():
    pts = [ [1,2,3,4,5,6,7,8,9],
            [2,3,5,7,3,6,8,9,1]]
    bspline = UniformBspline(pts)
    rg = bspline.getRegion()
    for i in range(rg[0],rg[1]):
        print '---------------'
        print i
        print bspline.evaluate(i+0.01)

if __name__ == '__main__':
    main()
         

        

     