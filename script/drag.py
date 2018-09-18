import numpy as np
import matplotlib.pyplot as plt
from uniform_bspline import UniformBspline

global pts
pts = [ [1,2,3,4,5,6,7,8,9],
        [2,3,5,7,3,6,8,9,1]]

global obss
obss = [[2,2],[3,3],[8,8],[9,5],[5,5],[7,3],[9,7]]
global obs

global move_id
move_id = None
global move_pos
move_pos = None

global fig, ax
fig, ax = plt.subplots()

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
    '''
    global pts
    bspline = UniformBspline(pts)
    rg = bspline.getRegion()
    u = rg[0]

    bx, by = [], []
    while u<rg[1]:
        val = bspline.evaluate(u)
        bx.append(val[0])
        by.append(val[1])
        u += 0.01

    cx, cy = [], []
    for i in range(rg[0]+1,rg[1]):
        val = bspline.evaluate(i)
        cx.append(val[0])
        cy.append(val[1])

    return bx, by, cx, cy


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
        x, y, cx, cy = drawBSpline()
        ax.clear()
        ax.axis([0,10,0,10])
        ax.plot(pts[0],pts[1],'ro')
        ax.plot(pts[0],pts[1],'g')
        ax.plot(x, y)
        ax.plot(cx, cy, 'yo')
        for i in range(len(obs)):
            ax.plot(obs[i][0],obs[i][1],'r')
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

    x, y, cx, cy = drawBSpline()
    ax.clear()
    ax.axis([0,10,0,10])
    ax.plot(pts[0],pts[1],'ro')
    ax.plot(pts[0],pts[1],'g')
    ax.plot(x, y)
    ax.plot(cx, cy, 'yo')
    for i in range(len(obs)):
        ax.plot(obs[i][0],obs[i][1],'r')
    plt.draw()


def main():
    
    global fig, ax, obs
    ax.axis([0,10,0,10])
    ax.plot(pts[0],pts[1],'ro')
    ax.plot(pts[0],pts[1],'g')
    x, y, cx, cy = drawBSpline()
    ax.plot(x, y)
    ax.plot(cx, cy, 'yo')
    drawObs()
    for i in range(len(obs)):
        ax.plot(obs[i][0],obs[i][1],'r')

    cid1 = fig.canvas.mpl_connect('button_press_event', onPress)
    cid2 = fig.canvas.mpl_connect('button_release_event', onRelease)
    cid3 = fig.canvas.mpl_connect('motion_notify_event', onMotion )

    plt.show()

if __name__ == '__main__':
    main()