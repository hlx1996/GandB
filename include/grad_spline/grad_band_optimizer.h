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

  public:
    GradBandOptimizer(Eigen::MatrixXd points, sdf_tools::SignedDistanceField* sdf, double res);
    ~GradBandOptimizer();

    // execute main operation here
    void optimize();

    // set signed distance field and resolution of it
    void setDistanceField(sdf_tools::SignedDistanceField* s, double res);

    // set the sequence points
    void setPoints(Eigen::MatrixXd points);
};
