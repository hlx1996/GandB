#include <ros/ros.h>
#include <visualization_msgs/Marker.h>

#include <Eigen/Eigen>

using namespace std;

// The state is a vector of {x, y, z, vx, vy, vz, ax, ay, az}
void stateTransit(Eigen::VectorXd init_state, Eigen::VectorXd& final_state, Eigen::Vector3d um, double tau)
{
    Eigen::Vector3d x0, v0, a0;
    x0 = init_state.head(3);
    v0 = init_state.block(3, 0, 3, 1);
    a0 = init_state.tail(3);

    // cout << "The initial state is:\n" << init_state.transpose() << endl;
    // cout << "um:" << um.transpose() << ", tau:" << tau << endl;

    // build the state transition matrix phi(tau)
    Eigen::MatrixXd phi(9, 9);
    phi.block(3, 0, 6, 6) = Eigen::MatrixXd::Zero(6, 6);
    phi.block(0, 0, 3, 3) = phi.block(6, 6, 3, 3) = Eigen::MatrixXd::Identity(3, 3);
    phi.block(0, 3, 3, 3) = phi.block(3, 6, 3, 3) = tau * Eigen::MatrixXd::Identity(3, 3);
    phi.block(0, 6, 3, 3) = Eigen::MatrixXd::Identity(3, 3) * 0.5 * tau * tau;

    // cout << "phi:\n" << phi << endl;

    // The build the integral terms
    Eigen::VectorXd integral(9);
    integral.head(3) = (1 / 3.0) * pow(tau, 3) * um;
    integral.segment(3, 3) = 0.5 * pow(tau, 2) * um;
    integral.tail(3) = tau * um;
    // cout << "integral:\n" << integral << endl;

    final_state = phi * init_state + integral;
    // cout << "The final state is:\n" << final_state.transpose() << endl << endl;
}

int main(int argc, char** argv)
{
    ros::init(argc, argv, "astar");
    ros::NodeHandle node;

    ros::Publisher primitive_pub = node.advertise<visualization_msgs::Marker>("primitive/path", 100);
    ros::Duration(0.5).sleep();

    cout << "test state transition" << endl;

    // given x0, v0, a0, we get the state by applying um for tau
    Eigen::VectorXd x0(9);
    x0 << 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0;

    Eigen::Vector3d um;
    um << 1.0, 1.0, 1.0;  // specify the input and duration
    double tau = 1.0;

    // calculate the final state
    Eigen::VectorXd x1;
    stateTransit(x0, x1, um, tau);

    // Then we can generate a random initial state, apply a range of um and tau and get all the motion primitive
    // We visualize all these primitive in rviz
    // At the beginning I don't care about the feasibility. So I just uniformly sample some um and tau

    double dum = 1.0;
    int idx = 0;
    for (int i = -1; i <= 1; ++i)
        for (int j = -1; j <= 1; ++j)
            for (int k = 0; k <= 0; ++k)
            {
                um << double(i) * dum, double(j) * dum, double(k) * dum;
                visualization_msgs::Marker p;
                p.header.frame_id = "world";
                p.header.stamp = ros::Time::now();
                p.id = idx;

                p.type = visualization_msgs::Marker::SPHERE_LIST;
                p.action = visualization_msgs::Marker::ADD;

                p.pose.orientation.x = 0.0;
                p.pose.orientation.y = 0.0;
                p.pose.orientation.z = 0.0;
                p.pose.orientation.w = 1.0;
                p.scale.x = p.scale.y = p.scale.z = 0.02;

                p.color.a = p.color.r = 1.0;
                p.color.g = p.color.b = 0.0;
                for (double l = 0.0; l < tau; l += tau / 100)
                {
                    stateTransit(x0, x1, um, l);
                    geometry_msgs::Point pt;
                    pt.x = x1.head(3)(0);
                    pt.y = x1.head(3)(1);
                    pt.z = x1.head(3)(2);
                    p.points.push_back(pt);
                }

                primitive_pub.publish(p);
                ros::Duration(0.001).sleep();
                ++idx;
            }

    // ros::Duration(10.0).sleep();

    return 0;
}