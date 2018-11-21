#include <geometry_msgs/PoseStamped.h>
#include <nav_msgs/Path.h>

#include <stdlib.h>
#include "display.h"

// sdf_tools
#include "sdf_tools/collision_map.hpp"
#include "sdf_tools/sdf.hpp"

#include "grad_spline/grad_band_optimizer.h"
#include "grad_spline/uniform_bspline.h"

using namespace std;

int main(int argc, char** argv)
{
    ros::init(argc, argv, "random");
    ros::NodeHandle node;

    // -------------------ros initialization---------------------
    setpoint_pub = node.advertise<visualization_msgs::Marker>("trajopt/setpoint", 10);
    traj_point_pub = node.advertise<visualization_msgs::Marker>("trajopt/traj_point", 10);
    traj_pub = node.advertise<nav_msgs::Path>("trajopt/init_traj", 5);
    ros::Publisher visualization_pub =
        node.advertise<visualization_msgs::Marker>("sdf_tools_tutorial_visualization", 1, true);
    va_pub = node.advertise<nav_msgs::Odometry>("trajopt/odom", 10);

    srand(ros::Time::now().toSec());
    ros::Duration(0.5).sleep();

    // parameter
    double obs_th, min_step, max_step;
    ros::param::get("/random/obs_th", obs_th);
    ros::param::get("/random/min_step", min_step);
    ros::param::get("/random/max_step", max_step);

    double alpha, beta, lamda1, lamda2, lamda3, dist0, st, scale;
    int num, point_opti_num, algorithm, solver, pow1, pow2, isloop;
    double max_vel, max_acc, interval;
    ros::param::get("/random/alpha", alpha);
    ros::param::get("/random/beta", beta);
    ros::param::get("/random/lamda1", lamda1);
    ros::param::get("/random/lamda2", lamda2);
    ros::param::get("/random/lamda3", lamda3);
    ros::param::get("/random/pow1", pow1);
    ros::param::get("/random/pow2", pow2);
    ros::param::get("/random/dist0", dist0);
    ros::param::get("/random/scale", scale);
    ros::param::get("/random/point_opti_num", point_opti_num);
    ros::param::get("/random/algorithm", algorithm);
    ros::param::get("/random/max_vel", max_vel);
    ros::param::get("/random/max_acc", max_acc);
    ros::param::get("/random/interval", interval);

    ros::param::get("/random/opti_num", num);
    ros::param::get("/random/sleep", st);
    ros::param::get("/random/solver", solver);
    ros::param::get("/random/loop", isloop);


    while (ros::ok())
    {
        // -------------------randomly generate some way points-------------------
        Eigen::Vector3d start, end;
        start(0) = start(1) = -4.0;
        end(0) = end(1) = 4.0;
        start(2) = end(2) = 2.0;

        vector<Eigen::Vector3d> way_points;
        way_points.push_back(start);

        double dist = (start - end).norm();
        for (int i = 0; dist > min_step;)
        {
            Eigen::Vector3d wp;
            wp(0) = way_points[i](0) + max_step * rand() / double(RAND_MAX);
            wp(1) = way_points[i](1) + max_step * rand() / double(RAND_MAX);
            wp(2) = 2.0;

            if (fabs(wp(0)) <= 4.0 && fabs(wp(1)) <= 4.0 && (wp - way_points[i]).norm() > min_step &&
                (wp - way_points[i]).norm() < max_step)
            {
                dist = (wp - end).norm();
                way_points.push_back(wp);
                ++i;
                // cout << "new way point:" << wp << endl;
            }
        }

        way_points.push_back(end);
        point_num = way_points.size();

        //  display  waypoints in rviz
        visualizeSetPoints(way_points);
        cout << "----------------------Way points created!----------------------" << endl;

        //---------------------create a map using sdf_tools-----------------------------

        // sdf collision map parameter
        const double resolution = 0.1;
        const double x_size = 20.0;
        const double z_size = 5.0;
        double y_size = 20.0;
        Eigen::Translation3d origin_translation(-10.0, -10.0, 0.0);
        Eigen::Quaterniond origin_rotation(1.0, 0.0, 0.0, 0.0);
        const Eigen::Isometry3d origin_transform = origin_translation * origin_rotation;
        const std ::string frame = "world";

        // create map
        sdf_tools ::COLLISION_CELL cell;
        cell.occupancy = 0.0;
        cell.component = 0;
        const sdf_tools ::COLLISION_CELL oob_cell = cell;
        sdf_tools ::CollisionMapGrid collision_map(origin_transform, frame, resolution, x_size, y_size, z_size,
                                                   oob_cell);

        // add some obstacle randomly
        sdf_tools::COLLISION_CELL obstacle_cell(1.0);

        int obs_num = 50;
        vector<Eigen::Vector3d> obstacles;
        cout << "----------------------Add obstacles to map!---------------" << endl;

        int fail_num = 0;
        for (int i = 0; i < obs_num;)
        {
            // randomly create a obstacle point
            Eigen::Vector3d pt;
            pt(0) = -3.5 + 7.0 * rand() / double(RAND_MAX);
            pt(1) = -3.5 + 7.0 * rand() / double(RAND_MAX);
            pt(2) = 2.0;

            // ensure that obstacle is far enough from start and end
            if ((pt - start).norm() < 1.0 || (pt - end).norm() < 1.0) continue;

            // ensure that obstacle is far enough from waypoint
            double dist_thresh = 0.4;
            double min_dist = 1000.0;
            for (int j = 0; j < way_points.size(); ++j)
            {
                double dist = (way_points[j] - pt).norm();
                if (dist < min_dist) min_dist = dist;
            }

            if (min_dist > dist_thresh)
            {
                // ensure that each obstacle is far enough from others
                if (i == 0)
                {
                    obstacles.push_back(pt);
                    ++i;
                }
                else
                {
                    min_dist = 1000.0;
                    dist_thresh = obs_th;
                    for (int j = 0; j < obstacles.size(); ++j)
                    {
                        double dist = (obstacles[j] - pt).norm();
                        if (dist < min_dist) min_dist = dist;
                    }

                    if (min_dist > dist_thresh)
                    {
                        obstacles.push_back(pt);
                        ++i;
                        fail_num = 0;
                    }
                    else
                    {
                        ++fail_num;
                    }
                }
            }
            if (fail_num > 10000)
            {
                break;
            }
        }

        cout << "----------------------Obstacles generated!----------------------" << endl;

        // add the generated obstacles into collision map
        const int th = 2;
        for (float z = 0; z < 3.5; z += resolution)
            for (int i = 0; i < obstacles.size(); ++i)
            {
                for (int m = -th; m <= th; m++)
                    for (int n = -th; n <= th; n++)
                    {
                        collision_map.Set(obstacles[i](0) + m * resolution, obstacles[i](1) + n * resolution, z,
                                          obstacle_cell);
                    }
            }

        // visualize the collision map
        std_msgs::ColorRGBA collision_color;
        collision_color.r = 0.0;
        collision_color.g = 0.0;
        collision_color.b = 1.0;
        collision_color.a = 0.8;

        std_msgs::ColorRGBA free_color, unknown_color;
        unknown_color.a = free_color.a = 0.0;

        visualization_msgs::Marker collision_map_marker =
            collision_map.ExportForDisplay(collision_color, free_color, unknown_color);
        collision_map_marker.ns = "collision_map";
        collision_map_marker.id = 1;

        visualization_pub.publish(collision_map_marker);

        // Build the signed distance field
        ros::Time t1 = ros::Time::now();

        float oob_value = INFINITY;
        std::pair<sdf_tools::SignedDistanceField, std::pair<double, double>> sdf_with_extrema =
            collision_map.ExtractSignedDistanceField(oob_value);
        sdf_tools::SignedDistanceField sdf = sdf_with_extrema.first;

        ros::Time t2 = ros::Time::now();

        cout << "Time for building distance field:" << (t2 - t1).toSec() << endl;
        cout << "----------------------Signed distance field build!----------------------" << endl;

        // ---------------------main optimization procedure---------------------
        // convert vector<Eigen::Vector3d> to Matrix Nx3
        Eigen::MatrixXd points = Eigen::MatrixXd::Zero(int(way_points.size()), 3);
        for (int i = 0; i < int(way_points.size()); ++i)
        {
            points(i, 0) = way_points.at(i)(0);
            points(i, 1) = way_points.at(i)(1);
            points(i, 2) = way_points.at(i)(2);
        }

        GradBandOptimizer optimizer(points, &sdf, resolution);
        optimizer.setParameter(alpha, beta, lamda1, lamda2, lamda3, pow1, pow2, dist0, scale, point_opti_num, algorithm,
                               max_vel, max_acc, interval);
        double time1 = 0.0;
        for (int i = 0; i < num; ++i)
        {
            ros::Duration(st).sleep();
            ros::Time t1 = ros::Time::now();
            if (solver == 1)
                optimizer.optimize();
            else if (solver == 2)
                optimizer.optimize2();
            else if (solver == 3)
                optimizer.optimize3();
            ros::Time t2 = ros::Time::now();
            time1 += t2.toSec() - t1.toSec();
            cout << "one time: " << t2 - t1 << " total time1:" << time1 << endl;
            Eigen::MatrixXd pts = optimizer.getPoints();
            visualizePoints(pts);

            // use the optimized points for bspline
            UniformBspline bspline = UniformBspline(pts, 4, interval);
            displayTrajectory(bspline);
        }

        if (!isloop) break;
        ros::Duration(2.0).sleep();
    }

    return 0;
}
