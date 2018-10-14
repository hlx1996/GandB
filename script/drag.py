import numpy as np
import matplotlib.pyplot as plt
from uniform_bspline import UniformBspline
from bspline import BSpline

global pts
pts = [ [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0],
        [2.0,3.0,5.0,7.0,3.0,6.0,8.0,9.0,1.0]]
# pts = [ [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0],
#         [5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0]]

global obss
obss = [[2,2],[3,3],[8,8],[9,5],[5,5],[7,3],[9,7]]
global obs

global move_id
move_id = None
global move_pos
move_pos = None

global fig, ax, ax0, ax1, ax2
fig = plt.figure()
ax = fig.add_subplot(221)
ax0 = fig.add_subplot(222)
ax1 = fig.add_subplot(223)
ax2 = fig.add_subplot(224)


def drawObs():
    '''
    A list of several obstales
    [[[],[]],[[],[]],[[],[]]...]

    '''
    global obss, obs
    obs = []
    for i in range(len(obss)):
        ob, obx, oby = [], [], []
        obx.append(obss[i][0])
        obx.append(obss[i][0])
        obx.append(obss[i][0]+1)
        obx.append(obss[i][0]+1)
        obx.append(obss[i][0])

        oby.append(obss[i][1])
        oby.append(obss[i][1]+1)
        oby.append(obss[i][1]+1)
        oby.append(obss[i][1])
        oby.append(obss[i][1])

        ob.append(obx)
        ob.append(oby)
        obs.append(ob)


def drawBSpline():
    '''
    Return the whole bspline and control point as bx, by, cx, cy
    Also return the velocity and acc and time as vx, vy, ax, ay, time
    '''
    global pts
    # bspline = UniformBspline(pts)
    bspline = BSpline(pts,4,False)
    rg = bspline.getRegion()
    d1 = bspline.getDerivative()
    d2 = d1.getDerivative()
    u = rg[0]

    bx, by = [], []
    vx, vy = [], []
    acx, acy = [], []
    time = []
    while u<rg[1]:
        # val = bspline.evaluate(u)
        val = bspline.evaluateByMat(u)
        bx.append(val[0])
        by.append(val[1])

        # val = bspline.evaluate(u,1)
        val = d1.evaluateByMat(u)
        vx.append(val[0])
        vy.append(val[1])


        # val = bspline.evaluate(u,2)
        val = d2.evaluateByMat(u)
        acx.append(val[0])
        acy.append(val[1])

        time.append(u)

        u += 0.01

    cx, cy = [], []
    for i in range(rg[0]+1,rg[1]):
        # val = bspline.evaluate(i,0)
        val = bspline.evaluateByMat(i)
        cx.append(val[0])
        cy.append(val[1])

    return bx, by, cx, cy, vx, vy, acx, acy, time


def onPress(event):
    '''
    If double click, add that point into the pts and redraw the bspline
    else move the clicked point 
    '''
    print('Press: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
          (event.button, event.x, event.y, event.xdata, event.ydata))

    if event.dblclick:
        print 'double click, add point into pts'
        clx, cly = event.xdata, event.ydata
        global pts
        for i in range(len(pts[0])):
            if pts[0][i] > clx:
                pts[0].insert(i,clx)
                pts[1].insert(i,cly)
                break

        'redraw'
        x, y, cx, cy, vx, vy, acx, acy, time = drawBSpline()
        global ax0, ax1, ax2
        ax.clear()
        ax0.clear()
        ax1.clear()
        ax2.clear()

        ax.axis([0,10,0,10])
        ax.plot(pts[0],pts[1],'ro')
        ax.plot(pts[0],pts[1],'g')
        ax.plot(x, y)
        ax.plot(cx, cy, 'yo')
        for i in range(len(obs)):
            ax.plot(obs[i][0],obs[i][1],'r')

        ax0.plot(time, x,'r')
        ax0.plot(time, y,'g')
        ax1.plot(time, vx,'r')
        ax1.plot(time, vy,'g')
        ax2.plot(time, acx,'r')
        ax2.plot(time, acy,'g')
        plt.draw()



    else:
        # get the closest point of pressed position, 
        x, y = event.xdata, event.ydata
        min_dst = 1000
        min_id = 0
        global pts
        for i in range(len(pts[0])):
            dst = np.sqrt((x-pts[0][i])**2+(y-pts[1][i])**2)
            if dst < min_dst:
                min_dst = dst
                min_id = i
        # print pts[0][min_id], pts[1][min_id]
        
        # The closest point is close enough?
        global move_id, move_pos
        if min_dst < 0.5:
            move_id = min_id
            move_pos = [x, y]
            print move_id, move_pos
        else:
            print 'no point press'
    
def onMotion(event):
    # print('Motion: x=%d, y=%d, xdata=%f, ydata=%f' %
    #       (event.x, event.y, event.xdata, event.ydata))
    global move_id, move_pos, pts
    if move_id is None: 
        # print 'None'
        return
    
    # calculate dx, dy to move_pos and update pts
    dx = event.xdata - move_pos[0]
    dy = event.ydata - move_pos[1]
    pts[0][move_id] += dx 
    pts[1][move_id] += dy 
    # print pts

    move_pos[0] = event.xdata
    move_pos[1] = event.ydata

    global ax
    # ax.clear()
    ax.axis([0,10,0,10])
    ax.plot(pts[0],pts[1],'ro')
    plt.draw()

        

def onRelease(event):
    # print('Release: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
    #       (event.button, event.x, event.y, event.xdata, event.ydata))
    global move_id, move_pos
    if move_id is None: 
        # print 'None'
        return
    move_id = None
    move_pos = None

    x, y, cx, cy, vx, vy, acx, acy, time = drawBSpline()

    global ax0, ax1, ax2
    ax.clear()
    ax0.clear()
    ax1.clear()
    ax2.clear()

    ax.axis([0,10,0,10])
    ax.plot(pts[0],pts[1],'ro')
    ax.plot(pts[0],pts[1],'g')
    ax.plot(x, y)
    ax.plot(cx, cy, 'yo')
    for i in range(len(obs)):
        ax.plot(obs[i][0],obs[i][1],'r')

    ax0.plot(time, x,'r')
    ax0.plot(time, y,'g')
    ax1.plot(time, vx,'r')
    ax1.plot(time, vy,'g')
    ax2.plot(time, acx,'r')
    ax2.plot(time, acy,'g')
    plt.draw()


def main():
    
    global fig, ax, obs, ax0, ax1, ax2
    ax.axis([0,10,0,10])
    ax.plot(pts[0],pts[1],'ro')
    ax.plot(pts[0],pts[1],'g')
    x, y, cx, cy, vx, vy, acx, acy, time = drawBSpline()
    ax.plot(x, y, 'r')
    ax.plot(cx, cy, 'yo')
    drawObs()
    for i in range(len(obs)):
        ax.plot(obs[i][0],obs[i][1],'r')

    ax0.plot(time, x,'r')
    ax0.plot(time, y,'g')
    ax1.plot(time, vx,'r')
    ax1.plot(time, vy,'g')
    ax2.plot(time, acx,'r')
    ax2.plot(time, acy,'g')

    cid1 = fig.canvas.mpl_connect('button_press_event', onPress)
    cid2 = fig.canvas.mpl_connect('button_release_event', onRelease)
    cid3 = fig.canvas.mpl_connect('motion_notify_event', onMotion )

    plt.show()

if __name__ == '__main__':
    main()
