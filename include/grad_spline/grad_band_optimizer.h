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
    double alpha, beta, lamda, dist0, alp, scale;
    int current_optimize_id, algorithm, point_opti_num;

    void getDistanceAndGradient(Eigen::Vector3d& pos, double& dist, Eigen::Vector3d& grad);

  public:
    GradBandOptimizer(Eigen::MatrixXd points, sdf_tools::SignedDistanceField* sdf, double res);
    ~GradBandOptimizer();

    // get optimized points
    Eigen::MatrixXd getPoints();

    // set algorithm parameters
    void setParameter(double alpha, double beta, double lamda, double dist0, double scale, int point_opti_num,
                      int algorithm);

    // execute main operation here
    // This function use hand-written optimization
    void optimize();
    // This function use NLopt optimization solver
    void optimize2();

    // set signed distance field and resolution of it
    void setDistanceField(sdf_tools::SignedDistanceField* s, double res);

    // set the sequence points
    void setPoints(Eigen::MatrixXd points);

    // For using NLopt solver, we need a func()
    static double costFunc(const std::vector<double>& x, std::vector<double>& grad, void* func_data);
};

#endif