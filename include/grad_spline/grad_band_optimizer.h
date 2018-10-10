#ifndef _GRAD_BAND_OPTIMIZER_H_
#define _GRAD_BAND_OPTIMIZER_H_

#include <Eigen/Eigen>

#include "sdf_tools/collision_map.hpp"
#include "sdf_tools/sdf.hpp"

// Gradient and elasitc band optimization

// Input: a signed distance field and a sequence of points
// Output: the optimized sequence of points
// The format of points: N x 3 matrix, each row is a point

class GradBandOptimizer
{
  private:
    /* data */
    sdf_tools::SignedDistanceField* sdf;
    double resolution;
    Eigen::MatrixXd points;
    double alpha, beta, lamda1, lamda2, lamda3, dist0, alp, scale;  // objective function parameters
    double max_vel, max_acc, interval;                              // constrains parameters
    int pow1, pow2;

    // used for adding derivative constrains
    // these seems not work... so try different way for adding costrains
    int current_optimize_id;  // id range from 0 to (points.rows-2)
    int current_axis;         // axis is 0, 1, 2 for x, y, z;
    int current_sign;         // sign is +1, -1 for lesser and larger
    int var_num;              // number of variables to be optimized

    int algorithm;
    int point_opti_num;

    void getDistanceAndGradient(Eigen::Vector3d& pos, double& dist, Eigen::Vector3d& grad);

  public:
    std::vector<double> min_var;
    double min_cost;

    GradBandOptimizer(Eigen::MatrixXd points, sdf_tools::SignedDistanceField* sdf, double res);
    ~GradBandOptimizer();

    // get optimized points
    Eigen::MatrixXd getPoints();

    // set algorithm parameters
    void setParameter(double alpha, double beta, double lamda1, double lamda2, double lamda3, int pow1, int pow2,
                      double dist0, double scale, int point_opti_num, int algorithm, double max_vel, double max_acc,
                      double interval);

    // execute main operation here
    // This function use hand-written optimization
    void optimize();
    // This function use NLopt optimization solver, and optimize each point in sequence
    void optimize2();
    // This function use NLopt optimization solver, and optimize all points in one time
    void optimize3();

    // set signed distance field and resolution of it
    void setDistanceField(sdf_tools::SignedDistanceField* s, double res);

    // set the sequence points
    void setPoints(Eigen::MatrixXd points);

    // For using NLopt solver, we need a func()
    static double costFunc2(const std::vector<double>& x, std::vector<double>& grad, void* func_data);

    static double costFunc3(const std::vector<double>& x, std::vector<double>& grad, void* func_data);
    static double costFunc4(const std::vector<double>& x, std::vector<double>& grad, void* func_data);
    static double costFunc5(const std::vector<double>& x, std::vector<double>& grad, void* func_data);
    static double velConstraint(const std::vector<double>& x, std::vector<double>& grad, void* data);
    static double accConstraint(const std::vector<double>& x, std::vector<double>& grad, void* data);
};

#endif