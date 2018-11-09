#ifndef _HYBRID_ASTAR_H
#define _HYBRID_ASTAR_H

#include <iostream>
#include <ros/ros.h>
#include <ros/console.h>
#include <Eigen/Eigen>
// #include "backward.hpp"
#include "grad_spline/data_type.h"
#include <string>
#include <map>

//#include <arc_utilities/voxel_grid.hpp>
#include <sdf_tools/collision_map.hpp>

using std::string;

class gridPathFinder
{
  private:
    Eigen::Vector3d gridIndex2coord(Eigen::Vector3i index);
    Eigen::Vector3i coord2gridIndex(Eigen::Vector3d pt);
    GridNodePtr pos2gridNodePtr(Eigen::Vector3d pos);

    double getDiagHeu(GridNodePtr node1, GridNodePtr node2);
    double getManhHeu(GridNodePtr node1, GridNodePtr node2);
    double getEuclHeu(GridNodePtr node1, GridNodePtr node2);
    double getHeu(GridNodePtr node1, GridNodePtr node2);
    bool shotHeu(GridNodePtr node1, GridNodePtr node2);

    void stateTransit1(Eigen::VectorXd init_state, Eigen::VectorXd& final_state, Eigen::Vector3d um, double tau);

    struct KinoState
    {
        double edge_cost;
        double heu;
        double optimal_time;
        double duration;
        Eigen::VectorXd input;
        Eigen::VectorXd state;
    };

    // diff, (edge_cost,heu_cost), (state, optimal_time)
    void getNeighbor(GridNodePtr current_node, GridNodePtr end_node, std::multimap<std::string, KinoState>& neighbors);

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
    gridPathFinder(Eigen::Vector3i GL_size, Eigen::Vector3i LOC_size)
    {
        // size of a big big global grid map
        GLX_SIZE = GL_size(0);
        GLY_SIZE = GL_size(1);
        GLZ_SIZE = GL_size(2);

        // size of local map, for local obs recording
        X_SIZE = LOC_size(0);
        Y_SIZE = LOC_size(1);
        Z_SIZE = LOC_size(2);
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

    // shot
    bool is_shot_succ = false;
    int N_max_shot = 10;
    int N_max = 10;
    int cnt_shot = 0;
    Eigen::MatrixXd coef_shot;
    double t_shot, dis_shot;

    bool has_path = false;

    GridNodePtr terminate_ptr;
};

#endif