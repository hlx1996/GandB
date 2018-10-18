// #include "grad_spline/a_star.h"
#include "grad_spline/hybrid_astar.h"
#include "display.h"
#include <ros/ros.h>

using namespace std;

int main(int argc, char** argv)
{
    ros::init(argc, argv, "astar");
    ros::NodeHandle node;

    ros::Publisher visualization_pub =
        node.advertise<visualization_msgs::Marker>("sdf_tools_tutorial_visualization", 1, true);

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
        new gridPathFinder(global_size, global_size);  // the global and local size are equal, which means the whole
                                                       // grid have occupancy info
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
        start(0) = start(1) = -4.0;
        end(0) = end(1) = 4.0;
        start(2) = end(2) = 2.0;

        int fail_num = 0;
        for (int i = 0; i < obs_num;)
        {
            // randomly create a obstacle point
            Eigen::Vector3d pt;
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
        path_finder->AstarSearch(start, end);
        vector<Eigen::Vector3d> path = path_finder->getPath();
        vector<GridNodePtr> path_nodes = path_finder->getPathNodes();
        vector<GridNodePtr> visited_nodes = path_finder->getVisitedNodes();
        path_finder->resetLocalMap();

        // displayPath(path, resolution);
        displayPathNodes(path_nodes, resolution);
        displayVisitedNodes(visited_nodes, resolution);

        delete collision_map;

        ros::Duration(0.1).sleep();

        break;
    }

    return 0;
}
