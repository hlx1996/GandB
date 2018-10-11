import rospy
from std_msgs.msg import String

from nav_msgs.msg import Odometry
import matplotlib.pyplot as plt

global vx
global vy
global vz 
global ax
global ay
global az
global t
global u
t = []
vx = []
vy = []
vz = []
ax = []
ay = []
az = []
u = 0.0

global data_ready
data_ready = False

global drawn
drawn = False

def callback(data):
    # rospy.loginfo(rospy.get_caller_id() + "I heard %s", data.data)
    global vx
    global vy
    global vz
    global ax
    global ay
    global az
    global t
    global u
    global drawn
    global data_ready

    if abs(data.pose.pose.position.x+1) < 1e-4 and abs(data.pose.pose.position.y+1) < 1e-4 and abs(data.pose.pose.position.z+1) < 1e-4:
        print "new data"
        t = []
        vx = []
        vy = []
        vz = []
        ax = []
        ay = []
        az = []
        u = 0.0
        data_ready = False
    elif abs(data.pose.pose.position.x-1) < 1e-4 and abs(data.pose.pose.position.y-1) < 1e-4 and abs(data.pose.pose.position.z-1) < 1e-4:
        data_ready = True
        print "ok to plot"
    else:
        vx.append(data.pose.pose.position.x)
        vy.append(data.pose.pose.position.y)
        vz.append(data.pose.pose.position.z)

        ax.append(data.twist.twist.linear.x)
        ay.append(data.twist.twist.linear.y)
        az.append(data.twist.twist.linear.z)

        t.append(u)
        u += 0.01
    
def listener():

    # In ROS, nodes are uniquely named. If two nodes with the same
    # node are launched, the previous one is kicked off. The
    # anonymous=True flag means that rospy will choose a unique
    # name for our 'listener' node so that multiple listeners can
    # run simultaneously.
    rospy.init_node('listener', anonymous=True)

    rospy.Subscriber("/trajopt/odom", Odometry, callback, queue_size=200)

    global data_ready
    global vx, vy, vz, ax, ay, az, t, drawn

    # use ion if want to redraw one by one
    plt.ion()
    plt.show()
    while not rospy.is_shutdown():
        print data_ready
        if data_ready:
            # draw
            print "draw"

            plt.clf()
            plt.figure(1, figsize=(12, 3))

            plt.subplot(121)
            plt.plot(t, vx, 'r')
            plt.plot(t, vy, 'g')
            plt.plot(t, vz, 'b')

            plt.subplot(122)
            plt.plot(t, ax, 'r')
            plt.plot(t, ay, 'g')
            plt.plot(t, az, 'b')
            plt.show()
            plt.pause(2.0)
            data_ready = False
            drawn = True
        else:
            rospy.sleep(0.5)

    # spin() simply keeps python from exiting until this node is stopped
    print "begin"
    rospy.spin()




if __name__ == '__main__':
    listener()
