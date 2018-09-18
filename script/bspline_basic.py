# This script is to test calculate the basic function of b-spline.
# Given the order p and n+1 control points, all the basic N_i_p are computed.
# These basic will be displayed can can be use later 

import numpy as np
import sympy as sp


def test():
    # x, y, z = sp.symbols("x y z")
    # a = x*(y+z)
    # sp.simplify(a)
    # print a 

    # a = np.zeros(5)
    # print a[1], a[2]

    a = range(5)
    b = np.zeros(5)
    for i in range(5):
        b[i] = a[4-i]
    print b

def calculateBasic():
    # basic n, p, m
    p = 4
    n = 50
    m = n + p + 1

    # compute knots vector
    u = []
    # for i in range(p+1):
    #     u.append(0)
    
    # for i in range(1,m-2*p):
    #     u.append(i)

    # for i in range(p+1):
    #     u.append(m-2*p)
    for i in range(m+1):
        u.append(i)

    print u

    # calculate the coefficient of Nip in each segment [u_k, u_k+1]
    v = sp.symbols("u")

    for k in range(p,m-p):
        coe = sp.Matrix(np.zeros(p+1))
        # print coe

        # 0 order Nip
        coe[0] = 1.0

        # 1-p order Nip
        for d in range(1,p+1):
            coe[d] = (v-u[k])/(u[k+d]-u[k])*coe[d-1]

            for i in range(1,d):
                coe[d-i] = coe[d-i-1]*(v-u[k-i])/(u[k+d-i]-u[k-i])+ coe[d-i]*(u[k+d-i+1]-v)/(u[k+d-i+1]-u[k-i+1])

            coe[0] = ((u[k+1]-v)/(u[k+1]-u[k-d+1]))*coe[0]

        coe = sp.expand(coe)

        # show the coefficient
        print 'The coefficient of [', u[k], ",", u[k+1], "] are:"
        print coe 








def main():
    test()
    calculateBasic()


if __name__ == '__main__':
    main() 