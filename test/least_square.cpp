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

    // define K(segment num), N(sample number), ts(segment time)
    double ts = 0.5;
    int K = 9;
    int N = 4;

    double tm = ts / N;

    // generate some samples, using sin(x), write bx, by, bz
    Eigen::VectorXd bx((N + 1) * (K + 1)), by((N + 1) * (K + 1)), bz((N + 1) * (K + 1));
    for (int i = 0; i <= K; ++i)
        for (int j = 0; j <= N; ++j)
        {
            double t = i * ts + j * tm;
            double x = t;
            double y = sin(2 * t);
            double z = 2.0;

            bx(i * (N + 1) + j) = x;
            by(i * (N + 1) + j) = y;
            bz(i * (N + 1) + j) = z;
        }
    cout << "bx:" << bx.transpose() << ", by:" << by.transpose() << ", bz:" << bz.transpose() << endl;

    // draw the origin curve
    double total_t = ts * (K + 1);
    vector<Eigen::Vector3d> path_origin;

    for (double t = 0.0; t <= total_t; t += 0.02)
    {
        Eigen::Vector3d pt;
        pt(0) = t;
        pt(1) = sin(2 * t);
        pt(2) = 2.0;
        path_origin.push_back(pt);
    }

    // draw the least square b-spline
    Eigen::MatrixXd samples(3, (K + 1) * (N + 1));
    samples.row(0) = bx;
    samples.row(1) = by;
    samples.row(2) = bz;
    Eigen::MatrixXd control_pts;
    getControlPointLeastSquare(K, N, ts, samples, control_pts);

    UniformBspline bspline(control_pts, 5, ts, false);

    vector<Eigen::Vector3d> path_bspline;
    for (double t = 0.0; t <= total_t; t += 0.02)
    {
        Eigen::Vector3d pt = bspline.evaluate(t);
        path_bspline.push_back(pt);
    }

    displayPathWithColor(path_origin, 0.05, 2, 0);
    displayPathWithColor(path_bspline, 0.05, 1, 1);
}