#include <iostream>
#include <Eigen/Dense>
#include <ros/ros.h>
#include "display.h"
#include "grad_spline/sdf_map.h"

int main(int argc, char** argv)
{
    ros::init(argc, argv, "sdf");
    ros::NodeHandle node;
    ros::Duration(1.0).sleep();

    ros::Publisher occ_pub = node.advertise<visualization_msgs::Marker>("sdf_map/occ", 1, true);
    ros::Publisher dist_pub = node.advertise<visualization_msgs::Marker>("sdf_map/dist", 1, true);

    ros::Time t1, t2;

    // create a map
    Eigen::Vector3d origin, map_size;
    origin << -5.0, -5.0, 0.0;
    map_size << 10.0, 10.0, 5.0;
    double resolution = 0.2;
    SDFMap sdf_map = SDFMap(origin, resolution, map_size);

    // set some pos's occ
    Eigen::Vector3d pos1;
    for (int x = -4; x <= 4; x += 2)
        for (int y = -4; y <= 4; y += 2)
            for (double z = 0.0; z < 5.0; z += 0.1)
            {
                if (x == 0 && y == 0) continue;
                pos1 << x, y, z;
                sdf_map.setOccupancy(pos1);
            }

    // visualize it
    visualization_msgs::Marker occ_marker;
    sdf_map.getOccupancyMarker(occ_marker, 0, Eigen::Vector4d(0, 0, 1.0, 0.8));
    occ_pub.publish(occ_marker);

    // update distance field
    // sdf_map.setUpdateRange();
    t1 = ros::Time::now();

    sdf_map.updateESDF3d();

    t2 = ros::Time::now();
    cout << "time in distance field:" << (t2 - t1).toSec() << endl;

    // visualize distance field
    vector<visualization_msgs::Marker> dis_markers;
    sdf_map.getESDFMarker(dis_markers, 0, Eigen::Vector3d(1.0, 0, 0));
    for (int i = 0; i < int(dis_markers.size()); ++i)
    {
        dist_pub.publish(dis_markers[i]);
        ros::Duration(0.1).sleep();
    }

    // get distance
    // cout << "distance at 000:" << sdf_map.getDistance(Eigen::Vector3d(0, 0, 0)) << endl;
    // cout << "distance at 110:" << sdf_map.getDistance(Eigen::Vector3d(1, 1, 0)) << endl;
    // cout << "distance at 220:" << sdf_map.getDistance(Eigen::Vector3d(2, 2, 0)) << endl;

    int input;
    cout << "wait for update" << endl;
    cin >> input;

    // add some new obstacle and update
    for (double z = 0.0; z < 5.0; z += 0.1)
    {
        pos1 << 0, 0, z;
        sdf_map.setOccupancy(pos1);
    }

    sdf_map.getOccupancyMarker(occ_marker, 0, Eigen::Vector4d(0, 0, 1.0, 0.8));
    occ_pub.publish(occ_marker);

    t1 = ros::Time::now();

    double max_dist = sdf_map.getMaxDistance();
    cout << "max dist:" << max_dist << endl;
    sdf_map.setUpdateRange(Eigen::Vector3d(0.0 - max_dist, 0.0 - max_dist, 0),
                           Eigen::Vector3d(0.0 + max_dist, 0.0 + max_dist, 5));
    sdf_map.updateESDF3d();

    t2 = ros::Time::now();
    cout << "time in distance field:" << (t2 - t1).toSec() << endl;

    // visualize distance field
    sdf_map.getESDFMarker(dis_markers, 0, Eigen::Vector3d(1.0, 0, 0));
    for (int i = 0; i < int(dis_markers.size()); ++i)
    {
        dist_pub.publish(dis_markers[i]);
        ros::Duration(0.1).sleep();
    }

    ros::Duration(1.0).sleep();

    return 0;
}