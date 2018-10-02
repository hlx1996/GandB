#ifndef _DISPLAY_H
#define _DISPLAY_H

#include <ros/ros.h>
#include <visualization_msgs/Marker.h>
#include <nav_msgs/Path.h>

#include <Eigen/Eigen>

#include <vector>
#include <stdlib.h>

#include "grad_spline/uniform_bspline.h"

using std::vector;

Eigen::VectorXd my_time;
int point_num;

ros::Publisher setpoint_pub;
ros::Publisher traj_pub;
ros::Publisher traj_point_pub;

// visualize initial waypoint
void visualizeSetPoints(vector<Eigen::Vector3d> points)
{
    // send them to rviz
    for (int i = 0; i < points.size(); ++i)
    {
        visualization_msgs::Marker p;
        p.header.frame_id = "world";
        p.header.stamp = ros::Time::now();
        p.id = i;

        p.type = visualization_msgs::Marker::SPHERE;
        p.action = visualization_msgs::Marker::ADD;

        p.pose.position.x = points[i][0];
        p.pose.position.y = points[i][1];
        p.pose.position.z = points[i][2];
        p.pose.orientation.w = 1;
        p.pose.orientation.x = 0;
        p.pose.orientation.y = 0;
        p.pose.orientation.z = 0;

        p.scale.x = p.scale.y = p.scale.z = 0.2;

        p.color.a = p.color.r = 1.0;
        p.color.g = p.color.b = 0.0;

        p.lifetime = ros::Duration(2000.0);

        setpoint_pub.publish(p);
        ros::Duration(0.001).sleep();

        // ROS_INFO_STREAM("publish set point");
    }
}

// visualize initial waypoint
void visualizePoints(Eigen::MatrixXd points)
{
    // send them to rviz
    srand(ros::Time::now().toSec());
    double cr = rand() / double(RAND_MAX);
    double cg = rand() / double(RAND_MAX);
    double cb = rand() / double(RAND_MAX);
    for (int i = 0; i < points.rows(); ++i)
    {
        visualization_msgs::Marker p;
        p.header.frame_id = "world";
        p.header.stamp = ros::Time::now();
        p.id = i;

        p.type = visualization_msgs::Marker::SPHERE;
        p.action = visualization_msgs::Marker::ADD;

        p.pose.position.x = points(i, 0);
        p.pose.position.y = points(i, 1);
        p.pose.position.z = points(i, 2);
        p.pose.orientation.w = 1;
        p.pose.orientation.x = 0;
        p.pose.orientation.y = 0;
        p.pose.orientation.z = 0;

        p.scale.x = p.scale.y = p.scale.z = 0.1;

        p.color.a = 1.0;
        p.color.r = cr;
        p.color.g = cg;
        p.color.b = cb;

        p.lifetime = ros::Duration(10.0);

        setpoint_pub.publish(p);
        ros::Duration(0.0001).sleep();

        // ROS_INFO_STREAM("publish set point");
    }
}

// use uniform bspline to draw trajectory
void displayTrajectory(UniformBspline bspline)
{
    nav_msgs::Path path;
    path.header.frame_id = "world";
    path.header.stamp = ros::Time::now();

    // publish these point
    double u1, u2;
    bspline.getRegion(u1, u2);
    double u = u1;

    while (u <= u2)
    {
        Eigen::Vector3d val = bspline.evaluate(u);
        geometry_msgs::PoseStamped pose;
        pose.header.frame_id = "world";

        pose.pose.position.x = val(0);
        pose.pose.position.y = val(1);
        pose.pose.position.z = val(2);

        pose.pose.orientation.w = 1;
        pose.pose.orientation.x = 0;
        pose.pose.orientation.y = 0;
        pose.pose.orientation.z = 0;

        path.poses.push_back(pose);

        u += 0.01;
    }

    traj_pub.publish(path);
}

#endif