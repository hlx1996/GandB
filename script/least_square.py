import numpy as np
import sympy as sp

def main():

    # control points
    p0 = sp.symbols("p0")
    p1 = sp.symbols("p1")
    p2 = sp.symbols("p2")
    p3 = sp.symbols("p3")
    p4 = sp.symbols("p4")
    p5 = sp.symbols("p5")

    # M6 matrix
    M6 = np.matrix([[1, 26, 66, 26, 1, 0], [-5, -50, 0, 50, 5, 0], [10, 20, -60, 20, 10, 0], [-10, 20, 0, -20, 10, 0], [5, -20, 30, -20, 5, 0], [-1, 5, -10, 10, -5, 1]])/(5.0*4*3*2) 
    print M6

    N = 4
    u_list = np.linspace(0, 1, N+1)
    print u_list

    pv = np.matrix([[p0], [p1], [p2], [p3], [p4], [p5]])

    for u in u_list:
        uv = np.matrix([1., u, u**2, u**3, u**4, u**5])
        print uv * M6 * pv

        u = 0.5
        uv = np.matrix([1., u, u**2, u**3, u**4, u**5])
        print uv * M6 * pv



    # s = sp.symbols("s")
    # # t = sp.symbols("T")
    # t = 1/s
    # dp = sp.symbols("dp")
    # v0 = sp.symbols("v0")
    # v1 = sp.symbols("v1")
    # wt = sp.symbols("wt")

    # # t = 1

    # delp = dp - v0 * t 
    # delv = v1 - v0 

    # # m = np.matrix([[t**5/120.0,t**4/24.0,t**3/6.0],[t**4/24.0,t**3/6.0,t**2/2.0],[t**3/6.0,t**2/2.0,t]])
    # m = np.matrix([[t**3/6.0,t**2/2.0],[t**2/2.0,t]])
    # n = np.matrix([[delp],[delv]])

    # a = sp.Matrix(m).inv()*sp.Matrix(n)

    # alpha = a[0]
    # beta = a[1]

    # print alpha
    # print beta

    # cost = 1/3.0* alpha**2 * t**3 + alpha * beta * t**2 + beta**2 * t + wt*t

    # cost = sp.collect(sp.simplify(cost), s)
    # print 
    # sp.pprint(cost)

    # diff = sp.collect(sp.simplify(sp.diff(cost, s)), s)
    # print 
    # sp.pprint(diff)



    

if __name__ == '__main__':
    main()

