import numpy as np
import sympy as sp

def main():
    s = sp.symbols("s")
    # t = sp.symbols("T")
    t = 1/s
    dp = sp.symbols("dp")
    v0 = sp.symbols("v0")
    v1 = sp.symbols("v1")
    wt = sp.symbols("wt")

    # t = 1

    delp = dp - v0 * t 
    delv = v1 - v0 

    # m = np.matrix([[t**5/120.0,t**4/24.0,t**3/6.0],[t**4/24.0,t**3/6.0,t**2/2.0],[t**3/6.0,t**2/2.0,t]])
    m = np.matrix([[t**3/6.0,t**2/2.0],[t**2/2.0,t]])
    n = np.matrix([[delp],[delv]])

    a = sp.Matrix(m).inv()*sp.Matrix(n)

    alpha = a[0]
    beta = a[1]

    print alpha
    print beta

    cost = 1/3.0* alpha**2 * t**3 + alpha * beta * t**2 + beta**2 * t + wt*t

    cost = sp.collect(sp.simplify(cost), s)
    print 
    sp.pprint(cost)

    diff = sp.collect(sp.simplify(sp.diff(cost, s)), s)
    print 
    sp.pprint(diff)


def cal():
    m1 = np.matrix([[1, 11, 11, 1]])
    m2 = np.matrix([[-1, 1, 0, 0, 0], [0, -1, 1, 0, 0], [0, 0, -1, 1, 0], [0, 0, 0, -1, 1]])

    print m1
    print m2

    print m1 * m2

    m1 = np.matrix([[1, 4, 1]])
    m2 = np.matrix([[1, -2, 1, 0, 0], [0, 1, -2, 1, 0], [0, 0, 1, -2, 1]])

    print m1
    print m2

    print m1 * m2
    

if __name__ == '__main__':
    # main()
    cal()

