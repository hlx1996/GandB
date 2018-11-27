#ifndef _HYBRID_ASTAR_H
#define _HYBRID_ASTAR_H

#include <ros/console.h>
#include <ros/ros.h>
#include <Eigen/Eigen>
#include <iostream>
// #include "backward.hpp"
#include <map>
#include <string>
#include "grad_spline/data_type.h"

//#include <arc_utilities/voxel_grid.hpp>
#include <sdf_tools/collision_map.hpp>

using std::string;

class gridPathFinder
{
  private:
    Eigen::MatrixXd phi;   // state transit matrix
    double max_tau = 0.3;  // transition time
    double max_vel = 3.0;
    double max_acc = 2.0;
    double max_jerk = 3.0;

    Eigen::Vector3d start_point, end_point;

    Eigen::Vector3d gridIndex2coord(Eigen::Vector3i index);
    Eigen::Vector3i coord2gridIndex(Eigen::Vector3d pt);
    GridNodePtr pos2gridNodePtr(Eigen::Vector3d pos);

    double getDiagHeu(GridNodePtr node1, GridNodePtr node2);
    double getManhHeu(GridNodePtr node1, GridNodePtr node2);
    double getEuclHeu(GridNodePtr node1, GridNodePtr node2);
    double getHeu(GridNodePtr node1, GridNodePtr node2);
    bool shotHeu(GridNodePtr node1, GridNodePtr node2);

    void stateTransit1(Eigen::VectorXd init_state, Eigen::VectorXd& final_state, Eigen::Vector3d um, double tau);
    void jerkInputStateTransit(Eigen::VectorXd init_state, Eigen::VectorXd& final_state, Eigen::Vector3d um,
                               double tau);

    struct KinoState
    {
        double edge_cost;
        double heu;
        double optimal_time;
        double duration;
        Eigen::VectorXd input;
        Eigen::VectorXd state;
    };

    class Neighbors
    {
      private:
        /* data */
        std::vector<std::pair<Eigen::Vector3i, KinoState>> neighbor_data;

      public:
        Neighbors(/* args */)
        {
        }
        ~Neighbors()
        {
        }

        void erase(Eigen::Vector3i id)
        {
            for (int i = 0; i < int(neighbor_data.size()); ++i)
            {
                if ((id - neighbor_data[i].first).norm() == 0)
                {
                    neighbor_data.erase(neighbor_data.begin() + i);
                    break;
                }
            }
        }

        void add(Eigen::Vector3i id, KinoState state)
        {
            neighbor_data.push_back(std::make_pair(id, state));
        }

        bool find(Eigen::Vector3i id, KinoState& state)
        {
            for (int i = 0; i < int(neighbor_data.size()); ++i)
            {
                if ((id - neighbor_data[i].first).norm() == 0)
                {
                    state = neighbor_data[i].second;
                    return true;
                }
            }

            return false;
        }

        std::vector<std::pair<Eigen::Vector3i, KinoState>> getData()
        {
            return this->neighbor_data;
        }
    };

    // diff, (edge_cost,heu_cost), (state, optimal_time)
    void getNeighbor(GridNodePtr current_node, GridNodePtr end_node, Neighbors& neighbors, double& vis_time);

    Eigen::Vector3i stringToIndex(std::string str);

    std::vector<GridNodePtr> retrievePath(GridNodePtr current);
    double getKinoDynamicHeu(GridNodePtr node1, GridNodePtr node2);
    double getKinoDynamicHeu(Eigen::VectorXd x1, Eigen::VectorXd x2, double& optimal_time);
    std::vector<double> quartic(double a, double b, double c, double d, double e);
    std::vector<double> cubic(double a, double b, double c, double d);

    double w_time;
    double resolution, inv_resolution;
    double gl_xl, gl_yl, gl_zl;
    double tie_breaker = 1.0 + 1.0 / 10000;

    std::vector<GridNodePtr> expandedNodes;
    std::vector<GridNodePtr> gridPath;

    int GLX_SIZE, GLY_SIZE, GLZ_SIZE;
    int X_SIZE, Y_SIZE, Z_SIZE;

    GridNodePtr*** GridNodeMap;
    std::multimap<double, GridNodePtr> openSet;

  public:
    gridPathFinder(Eigen::Vector3i GL_size, Eigen::Vector3i LOC_size, ros::NodeHandle node)
    {
        // size of a big big global grid map
        GLX_SIZE = GL_size(0);
        GLY_SIZE = GL_size(1);
        GLZ_SIZE = GL_size(2);

        // size of local map, for local obs recording
        X_SIZE = LOC_size(0);
        Y_SIZE = LOC_size(1);
        Z_SIZE = LOC_size(2);

        // state transit
        phi.resize(6, 6);
        phi.block(3, 0, 3, 3) = Eigen::MatrixXd::Zero(3, 3);
        phi.block(0, 0, 3, 3) = phi.block(3, 3, 3, 3) = Eigen::MatrixXd::Identity(3, 3);
        phi.block(0, 3, 3, 3) = Eigen::MatrixXd::Identity(3, 3);

        // path publisher
        path_pub = node.advertise<visualization_msgs::Marker>("hybridastar/path", 10);

        w_time = 5.0;
    };
    gridPathFinder(){};
    ~gridPathFinder(){};

    void initGridNodeMap(double _resolution, Eigen::Vector3d global_xyz_l);
    void linkLocalMap(sdf_tools::CollisionMapGrid* local_map, Eigen::Vector3d xyz_l);
    void AstarSearch(Eigen::Vector3d start_pt, Eigen::Vector3d start_vel, Eigen::Vector3d end_pt,
                     Eigen::Vector3d end_vel);

    void resetLocalMap();
    void resetPath();

    std::vector<Eigen::Vector3d> getPath();
    std::vector<GridNodePtr> getVisitedNodes();
    std::vector<GridNodePtr> getPathNodes();
    std::vector<Eigen::Vector3d> getKinoTraj(double resolution);
    Eigen::MatrixXd getSamples(double& ts, int& K, int N = 1, bool repeat = false);
    Eigen::MatrixXd getSamplesUniformLength(double ds, double& ts, int& K);

    // shot
    bool is_shot_succ = false;
    int N_max_shot = 10;
    int N_max = 10;
    int cnt_shot = 0;
    Eigen::MatrixXd coef_shot;
    double t_shot, dis_shot;

    bool has_path = false;

    GridNodePtr terminate_ptr;

    ros::Publisher path_pub;

    // color 1234, rgby, resolution = 0.05
    void displayPathWithColor(std::vector<Eigen::Vector3d> path, double resolution, int color, int id)
    {
        visualization_msgs::Marker mk;
        mk.header.frame_id = "world";
        mk.header.stamp = ros::Time::now();
        mk.type = visualization_msgs::Marker::SPHERE_LIST;
        mk.action = visualization_msgs::Marker::DELETE;
        mk.id = id;

        path_pub.publish(mk);

        mk.action = visualization_msgs::Marker::ADD;
        mk.pose.orientation.x = 0.0;
        mk.pose.orientation.y = 0.0;
        mk.pose.orientation.z = 0.0;
        mk.pose.orientation.w = 1.0;
        mk.color.a = 1.0;

        if (color == 1)
        {
            mk.color.r = 1.0;
            mk.color.g = 0.0;
            mk.color.b = 0.0;
        }
        else if (color == 2)
        {
            mk.color.g = 1.0;
            mk.color.r = 0.0;
            mk.color.b = 0.0;
        }
        else if (color == 3)
        {
            mk.color.b = 1.0;
            mk.color.g = 0.0;
            mk.color.r = 0.0;
        }
        else if (color == 4)
        {
            mk.color.g = 1.0;
            mk.color.r = 1.0;
            mk.color.b = 0.0;
        }

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
};

#endif