// #include "grad_spline/a_star.h"
#include "grad_spline/hybrid_astar.h"
#include "grad_spline/uniform_bspline.h"
#include "display.h"
#include <ros/ros.h>

using namespace std;

int main(int argc, char** argv)
{
    ros::init(argc, argv, "astar");
    ros::NodeHandle node;

    ros::Publisher visualization_pub =
        node.advertise<visualization_msgs::Marker>("sdf_tools_tutorial_visualization", 1, true);
    va_pub = node.advertise<nav_msgs::Odometry>("trajopt/odom", 10);

    path_pub = node.advertise<visualization_msgs::Marker>("astar/path", 10);
    visited_pub = node.advertise<visualization_msgs::Marker>("astar/visited", 10);

    ros::Duration(1.0).sleep();

    srand(ros::Time::now().toSec());

    // create a a_star path finder
    double resolution = 0.2;
    int map_size = 10;
    int index_size = map_size / resolution;
    Eigen::Vector3i global_size(index_size, index_size, index_size);
    gridPathFinder* path_finder =
        new gridPathFinder(global_size, global_size, node);  // the global and local size are equal, which means the
                                                             // whole grid have occupancy info
    Eigen::Vector3d origin;
    origin << -map_size / 2.0, -map_size / 2.0, 0.0;
    path_finder->initGridNodeMap(resolution, origin);

    while (ros::ok())
    {
        // create a collision map for the path finder
        Eigen::Translation3d origin_translation(-map_size / 2.0, -map_size / 2.0, 0.0);
        Eigen::Quaterniond origin_rotation(1.0, 0.0, 0.0, 0.0);
        const Eigen::Isometry3d origin_transform = origin_translation * origin_rotation;
        sdf_tools::COLLISION_CELL free_cell(0.0);
        sdf_tools::CollisionMapGrid* collision_map = new sdf_tools::CollisionMapGrid(
            origin_transform, "world", resolution, map_size, map_size, map_size, free_cell);

        // the initial collision map is free, so we should add some random obstacle
        sdf_tools::COLLISION_CELL obstacle_cell(1.0);

        int obs_num = 50;
        vector<Eigen::Vector3d> obstacles;
        Eigen::Vector3d start, end;
        start(0) = start(1) = -4.0 + 1e-3;
        end(0) = 4.0 + 1e-3;
        end(1) = 4.0 + 1e-3;
        start(2) = end(2) = 2.0 + 1e-3;

        // add a obstacle near start
        Eigen::Vector3d pt;
        pt << -3.5, -3.5, 2.0;
        obstacles.push_back(pt);

        int fail_num = 0;
        for (int i = 0; i < obs_num;)
        {
            // randomly create a obstacle point
            pt(0) = -3.5 + 7.0 * rand() / double(RAND_MAX);
            pt(1) = -3.5 + 7.0 * rand() / double(RAND_MAX);
            pt(2) = 2.0;

            // ensure that obstacle is far enough from start and end
            if ((pt - start).norm() < 2.0 || (pt - end).norm() < 2.0) continue;

            // ensure that each obstacle is far enough from others
            if (i == 0)
            {
                obstacles.push_back(pt);
                ++i;
            }
            else
            {
                double min_dist = 1000.0;
                static double dist_thresh = 1.6;
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
                    ++fail_num;
            }
            if (fail_num > 10000) break;
        }

        // add the generated obstacles into collision map
        const int th = 1;
        for (float z = 0; z < 3.5; z += resolution)
            for (int i = 0; i < obstacles.size(); ++i)
                for (int m = -th; m <= th; m++)
                    for (int n = -th; n <= th; n++)
                    {
                        collision_map->Set(obstacles[i](0) + m * resolution, obstacles[i](1) + n * resolution, z,
                                           obstacle_cell);
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
            collision_map->ExportForDisplay(collision_color, free_color, unknown_color);
        collision_map_marker.ns = "collision_map";
        collision_map_marker.id = 1;
        visualization_pub.publish(collision_map_marker);

        // finally we finish the map... add to path finder
        path_finder->linkLocalMap(collision_map, origin);

        // now we can do the path searching, get the path and visited node
        Eigen::Vector3d sv, ev;
        sv << 0.0, 1.5, 0.0;
        ev << 1.5, 2.0, 0.0;
        path_finder->AstarSearch(start, sv, end, ev);
        // vector<Eigen::Vector3d> path = path_finder->getPath();
        vector<GridNodePtr> path_nodes = path_finder->getPathNodes();
        vector<GridNodePtr> visited_nodes = path_finder->getVisitedNodes();
        vector<Eigen::Vector3d> kino_path = path_finder->getKinoTraj(0.01);

        displayPathWithColor(kino_path, resolution * 0.25, 1, 1);
        displayVisitedNodes(visited_nodes, resolution);

        // display the start and end points
        vector<Eigen::Vector3d> start_end;
        start_end.push_back(start);
        start_end.push_back(end);
        // displayPathWithColor(start_end, 0.15, 2, 2);

        delete collision_map;

        // convert the path into b-spline using least square, N = 4,ts = 0.5
        ros::Time t1 = ros::Time::now();

        int K;
        double ts = 0.10;
        // get sample
        Eigen::MatrixXd samples = path_finder->getSamples(ts, K);
        // Eigen::MatrixXd samples = path_finder->getSamplesUniformLength(ds, ts, K);

        ros::Time t2 = ros::Time::now();
        cout << "time in get samples:" << (t2 - t1).toSec() << endl;

        vector<Eigen::Vector3d> v_sample;
        for (int i = 0; i < int(samples.cols()); ++i)
        {
            Eigen::Vector3d sample = samples.block(0, i, 3, 1);
            v_sample.push_back(sample);
        }
        displayPathWithColor(v_sample, 0.1, 4, 3);

        // convert to b-spline
        t1 = ros::Time::now();

        Eigen::MatrixXd control_pts;
        Eigen::MatrixXd extend_samples;
        extend_samples.resize(3, K + 6);
        extend_samples.block(0, 0, 3, K + 2) = samples;
        extend_samples.col(K + 2) = sv;
        extend_samples.col(K + 3) = ev;
        extend_samples.col(K + 4) = Eigen::Vector3d::Zero();
        extend_samples.col(K + 5) = Eigen::Vector3d::Zero();
        getControlPointEqu(extend_samples, ts, control_pts);
        // getControlPointLeastSquare(K, N, ts, samples.block(1, 0, 3, (N + 1) * (K + 1)), control_pts);

        t2 = ros::Time::now();
        cout << "time in Ax=b:" << (t2 - t1).toSec() << endl;

        // find the feasible ts
        t1 = ros::Time::now();
        while (true)
        {
            // get control point of vel acc
            UniformBspline bspline(control_pts, 5, ts, false);
            UniformBspline vt = bspline.getDerivative();
            UniformBspline at = vt.getDerivative();

            Eigen::MatrixXd vc, ac;
            vc = vt.getControlPoint();
            ac = at.getControlPoint();

            bool feasible = true;

            // check vel
            for (int i = 0; i < int(vc.rows()); ++i)
            {
                Eigen::Vector3d vci = vc.row(i);
                double v_max = max(fabs(vci(0)), fabs(vci(1)));
                v_max = max(fabs(vci(2)), v_max);

                if (v_max > 5.0)
                {
                    ts *= 1.1;
                    feasible = false;
                    break;
                }
            }
            if (!feasible) continue;

            // check acc
            for (int i = 0; i < int(ac.rows()); ++i)
            {
                Eigen::Vector3d aci = ac.row(i);
                double a_max = max(fabs(aci(0)), fabs(aci(1)));
                a_max = max(fabs(aci(2)), a_max);

                if (a_max > 5.0)
                {
                    ts *= 1.1;
                    feasible = false;
                    break;
                }
            }
            if (!feasible) continue;

            break;
        }
        cout << "new ts:" << ts << ", total time:" << (K + 1) * ts << endl;

        t2 = ros::Time::now();

        // draw the bspline
        UniformBspline bspline(control_pts, 5, ts, false);
        vector<Eigen::Vector3d> path_bspline;
        double total_t = ts * (K + 1);
        for (double t = 0.0; t <= total_t; t += 0.02)
        {
            Eigen::Vector3d pt = bspline.evaluate(t);
            path_bspline.push_back(pt);
        }
        displayPathWithColor(path_bspline, 0.05, 3, 3);

        // draw the control point
        vector<Eigen::Vector3d> ctp;
        for (int i = 0; i < int(control_pts.rows()); ++i)
        {
            Eigen::Vector3d pt = control_pts.row(i).transpose();
            ctp.push_back(pt);
        }
        displayPathWithColor(ctp, 0.15, 4, 4);
        cout << "control pt num:" << ctp.size() << endl;

        // draw the first and last segment control point
        ctp.erase(ctp.begin() + 6, ctp.end() - 6);
        displayPathWithColor(ctp, 0.25, 5, 1);

        // draw the vel and acc of b-spline
        UniformBspline velspline = bspline.getDerivative();
        UniformBspline accspline = velspline.getDerivative();
        vector<Eigen::Vector3d> vels, accs;
        for (double t = 0.0; t <= total_t; t += 0.01)
        {
            Eigen::Vector3d vt = velspline.evaluate(t);
            Eigen::Vector3d at = accspline.evaluate(t);
            vels.push_back(vt);
            accs.push_back(at);

            if (fabs(vt(0)) > 3.0 || fabs(vt(1)) > 3.0 || fabs(vt(2)) > 3.0) cout << "vel outside" << endl;
            if (fabs(at(0)) > 3.0 || fabs(at(1)) > 3.0 || fabs(at(2)) > 3.0) cout << "acc outside" << endl;
        }
        displayPathWithColor(vels, 0.1, 2, 5);
        displayPathWithColor(accs, 0.1, 3, 6);

        // clear path finder
        path_finder->resetLocalMap();

        ros::Duration(0.1).sleep();

        break;
    }

    return 0;
}
