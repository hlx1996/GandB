#include <ros/ros.h>
#include <visualization_msgs/Marker.h>

#include <Eigen/Eigen>

#include <string>
#include <map>
#include <sstream>
#include "display.h"

using namespace std;

ros::Publisher primitive_pub;
ros::Publisher neighbor_pub;

int idx = 0;

// The state is a vector of {x, y, z, vx, vy, vz, ax, ay, az}
void stateTransit2(Eigen::VectorXd init_state, Eigen::VectorXd& final_state, Eigen::Vector3d um, double tau)
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

// The state is a vector of {x, y, z, vx, vy, vz}
void stateTransit1(Eigen::VectorXd init_state, Eigen::VectorXd& final_state, Eigen::Vector3d um, double tau)
{
    Eigen::Vector3d x0, v0;
    x0 = init_state.head(3);
    v0 = init_state.tail(3);

    // cout << "The initial state is:\n" << init_state.transpose() << endl;
    // cout << "um:" << um.transpose() << ", tau:" << tau << endl;

    // build the state transition matrix phi(tau)
    Eigen::MatrixXd phi(6, 6);
    phi.block(3, 0, 3, 3) = Eigen::MatrixXd::Zero(3, 3);
    phi.block(0, 0, 3, 3) = phi.block(3, 3, 3, 3) = Eigen::MatrixXd::Identity(3, 3);
    phi.block(0, 3, 3, 3) = tau * Eigen::MatrixXd::Identity(3, 3);

    // cout << "phi:\n" << phi << endl;

    // The build the integral terms
    Eigen::VectorXd integral(6);
    integral.head(3) = 0.5 * pow(tau, 2) * um;
    integral.tail(3) = tau * um;
    // cout << "integral:\n" << integral << endl;

    final_state = phi * init_state + integral;
    // cout << "The final state is:\n" << final_state.transpose() << endl << endl;
}

void jerkInputStateTransit(Eigen::VectorXd init_state, Eigen::VectorXd& final_state, Eigen::Vector3d um, double tau)
{
    Eigen::Vector3d x0, v0, a0;
    x0 = init_state.head(3);
    v0 = init_state.segment(3, 3);
    a0 = init_state.tail(3);

    // add tau in transit matrix
    Eigen::MatrixXd phi(9, 9);
    phi.block(3, 0, 3, 3) = phi.block(6, 0, 3, 3) = phi.block(6, 3, 3, 3) = Eigen::MatrixXd::Zero(3, 3);
    phi.block(0, 0, 3, 3) = phi.block(3, 3, 3, 3) = phi.block(6, 6, 3, 3) = Eigen::MatrixXd::Identity(3, 3);
    phi.block(0, 3, 3, 3) = phi.block(3, 6, 3, 3) = Eigen::MatrixXd::Identity(3, 3);
    phi.block(0, 6, 3, 3) = 0.5 * Eigen::MatrixXd::Identity(3, 3);
    for (int i = 0; i < 6; ++i) phi(i, i + 3) = tau;
    for (int i = 0; i < 3; ++i) phi(i, i + 6) = 0.5 * tau * tau;
    // std::cout << "phi:\n" << phi << std::endl;

    // integral
    Eigen::VectorXd integral(9);
    integral.head(3) = (1 / 6.0) * pow(tau, 3) * um;
    integral.segment(3, 3) = 0.5 * pow(tau, 2) * um;
    integral.tail(3) = tau * um;
    // std::cout << "integral:\n" << integral << std::endl;

    // cout << "init:" << init_state << endl;
    final_state = phi * init_state + integral;
    // cout << "final: " << final_state.transpose() << endl;
}

int main(int argc, char** argv)
{
    ros::init(argc, argv, "astar");
    ros::NodeHandle node;

    path_pub = node.advertise<visualization_msgs::Marker>("astar/path", 10);

    ros::Duration(0.5).sleep();

    double max_jerk;
    double max_tau;
    cout << "input jerk ,tau" << endl;
    cin >> max_jerk >> max_tau;
    cout << max_jerk << ", " << max_tau << endl;

    Eigen::Vector3d um;
    double res = 1 / 2.0, time_res = 1 / 100.0;

    Eigen::VectorXd state(9);
    state << 0, 0, 0, 2, 2, 0, 0, 0, 0;
    int id = 0;

    for (double jx = -max_jerk; jx <= max_jerk + 1e-3; jx += max_jerk * res)
        for (double jy = -max_jerk; jy <= max_jerk + 1e-3; jy += max_jerk * res)
            for (double jz = -1e-3; jz <= 0.0 + 1e-3; jz += max_jerk * res)
            {
                um << jx, jy, jz;
                cout << "um:" << um.transpose() << endl;

                vector<Eigen::Vector3d> path;
                Eigen::VectorXd x1;
                for (double tau = max_tau * time_res; tau <= max_tau + 1e-3; tau += max_tau * time_res)
                {
                    jerkInputStateTransit(state, x1, um, tau);
                    path.push_back(x1.head(3));
                }

                displayPathWithColor(path, 0.02, 1, id);
                ++id;

                jerkInputStateTransit(state, x1, um, max_tau);
                cout << "state:" << x1.transpose() << endl;
            }

    return 0;
}