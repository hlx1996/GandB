#include "grad_spline/grad_band_optimizer.h"
using namespace std;

GradBandOptimizer::GradBandOptimizer(Eigen::MatrixXd points, sdf_tools::SignedDistanceField* sdf, double res)
{
    this->sdf = sdf;
    this->resolution = res;
    this->points = points;

    cout << "points:\n" << this->points << endl;
}

GradBandOptimizer::~GradBandOptimizer()
{
}

void GradBandOptimizer::optimize()
{
    cout << "begin optimization" << endl;

    // The first and last points remain unchanged
    // Other intermediate points have three forces: Two from neighbor points and one from gradient
    // F1 = Q_i+1 - Q_i, F2 = Q_i-1 - Q_i, F3 = -(p-p0)*Grad or 0
    // F = F1 + F2 + F3, delta_q = alpha * F
    // Q += delta_q
    // alpha *= beta


}

void GradBandOptimizer::setDistanceField(sdf_tools::SignedDistanceField* s, double res)
{
    this->sdf = sdf;
    this->resolution = res;
}

void GradBandOptimizer::setPoints(Eigen::MatrixXd points)
{
    this->points = points;
}