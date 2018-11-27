#include <iostream>
#include <Eigen/Dense>
#include <ros/ros.h>
#include "grad_spline/uniform_bspline.h"
#include "display.h"

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{
    ros::init(argc, argv, "least_square");
    ros::NodeHandle node;
    path_pub = node.advertise<visualization_msgs::Marker>("astar/path", 10);
    ros::Duration(1.0).sleep();

    int K = 50;  // K+1 segments, K+2 end points, K+6 control pts(constrain)
    double T = 8.0;
    double ts = T / (K + 1);

    // generate a sin curve in [0, 4]
    vector<Eigen::Vector3d> curve;
    double pi = 3.141592;
    for (double t = 0.0; t <= T; t += 0.01)
    {
        Eigen::Vector3d pt;
        pt << t - 4, sin(0.5 * pi * (t - 4)), 2.0;
        curve.push_back(pt);
    }
    displayPathWithColor(curve, 0.05, 1, 0);

    // get K+2 pos samples
    Eigen::MatrixXd samples(3, K + 6);
    int k = 0;
    for (double t = 0.0; t <= T + 1e-3; t += ts)
    {
        // cout << "k:" << k << endl;
        samples(0, k) = t - 4;
        samples(1, k) = sin(0.5 * pi * (t - 4));
        samples(2, k) = 2.05;
        ++k;
    }
    // cout << "final k:" << k << endl;

    // get vel samples
    samples(0, K + 2) = 1;
    samples(1, K + 2) = 0.5 * pi * cos(0.5 * pi * (0 - 4));
    samples(2, K + 2) = 0.0;

    samples(0, K + 3) = 1;
    samples(1, K + 3) = 0.5 * pi * cos(0.5 * pi * (T - 4));
    samples(2, K + 3) = 0.0;

    // get acc samples
    samples(0, K + 4) = 0.0;
    // samples by(K + 4) = -0.25 * pi * pi * sin(0.5 * pi * (0 - 4));
    samples(1, K + 4) = 0.0;
    samples(2, K + 4) = 0.0;

    samples(0, K + 5) = 0.0;
    samples(1, K + 5) = 0.0;
    samples(2, K + 5) = 0.0;

    // solve the equation
    // cout << "solve:" << endl;
    ros::Time t1 = ros::Time::now();

    Eigen::MatrixXd control_pts;
    getControlPointEqu(samples, ts, control_pts);

    ros::Time t2 = ros::Time::now();
    cout << "time:" << (t2 - t1).toSec() << endl;

    // draw bspline

    UniformBspline bspline(control_pts, 5, ts, false);

    vector<Eigen::Vector3d> path_bspline;
    for (double t = 0.0; t <= T; t += 0.01)
    {
        Eigen::Vector3d pt = bspline.evaluate(t);
        path_bspline.push_back(pt);
    }
    displayPathWithColor(path_bspline, 0.05, 2, 1);

    // draw control points
    vector<Eigen::Vector3d> ctps;
    for (int i = 0; i < control_pts.rows(); ++i)
    {
        Eigen::Vector3d pt = control_pts.row(i);
        ctps.push_back(pt);
    }
    displayPathWithColor(ctps, 0.1, 4, 2);

    ros::Duration(1.0).sleep();

    return 0;
}