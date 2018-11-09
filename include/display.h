#ifndef _DISPLAY_H
#define _DISPLAY_H

#include <ros/ros.h>
#include <visualization_msgs/Marker.h>
#include <nav_msgs/Path.h>
#include <nav_msgs/Odometry.h>

#include <Eigen/Eigen>

#include <vector>
#include <stdlib.h>

#include "grad_spline/uniform_bspline.h"
#include "grad_spline/data_type.h"

using std::vector;

Eigen::VectorXd my_time;
int point_num;

ros::Publisher setpoint_pub;
ros::Publisher traj_pub;
ros::Publisher traj_point_pub;
ros::Publisher va_pub;
ros::Publisher path_pub;
ros::Publisher visited_pub;

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

        p.lifetime = ros::Duration(2000.0);

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
    UniformBspline vel = bspline.getDerivative();
    UniformBspline acc = vel.getDerivative();

    ros::Time start = ros::Time::now();

    // start msg
    nav_msgs::Odometry odom;
    odom.header.frame_id = "world";
    odom.header.stamp = ros::Time::now();
    odom.pose.pose.position.x = -1.0;
    odom.pose.pose.position.y = -1.0;
    odom.pose.pose.position.z = -1.0;
    va_pub.publish(odom);

    while (u <= u2)
    {
        // publish the whole trajectory
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

        // publish the velocity and acceleration
        Eigen::Vector3d v = vel.evaluate(u);
        Eigen::Vector3d a = acc.evaluate(u);

        nav_msgs::Odometry odom;
        odom.header.frame_id = "world";
        odom.header.stamp = ros::Time(u);
        // odom.header.stamp = ros::Time::now();
        odom.pose.pose.position.x = v(0);
        odom.pose.pose.position.y = v(1);
        odom.pose.pose.position.z = v(2);
        odom.twist.twist.linear.x = a(0);
        odom.twist.twist.linear.y = a(1);
        odom.twist.twist.linear.z = a(2);
        va_pub.publish(odom);
        ros::Duration(0.001).sleep();

        u += 0.01;
    }

    odom.header.frame_id = "world";
    odom.header.stamp = ros::Time::now();
    odom.pose.pose.position.x = 1.0;
    odom.pose.pose.position.y = 1.0;
    odom.pose.pose.position.z = 1.0;
    va_pub.publish(odom);

    traj_pub.publish(path);
}

void displayPath(vector<Eigen::Vector3d> path, double resolution)
{
    visualization_msgs::Marker mk;
    mk.header.frame_id = "world";
    mk.header.stamp = ros::Time::now();
    mk.type = visualization_msgs::Marker::SPHERE_LIST;
    mk.action = visualization_msgs::Marker::ADD;
    mk.id = 0;

    mk.pose.orientation.x = 0.0;
    mk.pose.orientation.y = 0.0;
    mk.pose.orientation.z = 0.0;
    mk.pose.orientation.w = 1.0;
    mk.color.a = 1.0;
    mk.color.r = 1.0;
    mk.color.g = 0.0;
    mk.color.b = 0.0;

    mk.scale.x = resolution;
    mk.scale.y = resolution;
    mk.scale.z = resolution;

    geometry_msgs::Point pt;
    for (int i = 0; i < int(path.size()); i++)
    {
        pt.x = path[i](0);
        pt.y = path[i](1);
        pt.z = path[i](2);
        mk.points.push_back(pt);
    }
    path_pub.publish(mk);
    ros::Duration(0.001).sleep();
}

void displayVisitedNodes(vector<GridNodePtr> nodes, double resolution)
{
    visualization_msgs::Marker mk;
    mk.header.frame_id = "world";
    mk.header.stamp = ros::Time::now();
    mk.type = visualization_msgs::Marker::CUBE_LIST;
    mk.action = visualization_msgs::Marker::ADD;
    mk.id = 0;

    mk.pose.orientation.x = 0.0;
    mk.pose.orientation.y = 0.0;
    mk.pose.orientation.z = 0.0;
    mk.pose.orientation.w = 1.0;
    mk.color.a = 0.15;
    mk.color.r = 0.0;
    mk.color.g = 1.0;
    mk.color.b = 0.0;

    mk.scale.x = resolution;
    mk.scale.y = resolution;
    mk.scale.z = resolution;

    geometry_msgs::Point pt;
    for (int i = 0; i < int(nodes.size()); i++)
    {
        Eigen::Vector3d coord = nodes[i]->coord;
        pt.x = coord(0);
        pt.y = coord(1);
        pt.z = coord(2);

        mk.points.push_back(pt);
    }

    visited_pub.publish(mk);
}

#endif