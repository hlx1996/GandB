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

Eigen::MatrixXd GradBandOptimizer::getPoints()
{
    return this->points;
}

void GradBandOptimizer::setParameter(double alpha, double beta, double lamda, double dist0)
{
    this->alpha = alpha;
    this->alp = alpha;
    this->beta = beta;
    this->lamda = lamda;
    this->dist0 = dist0;
}

void GradBandOptimizer::getDistanceAndGradient(Eigen::Vector3d& pos, double& dist, Eigen::Vector3d& grad)
{
    std::vector<double> gradient = this->sdf->GetGradient(pos(0), pos(1), pos(2), true);
    grad(0) = gradient[0];
    grad(1) = gradient[1];
    grad(2) = gradient[2];
    std::pair<float, bool> location = this->sdf->GetSafe(pos(0), pos(1), pos(2));
    dist = location.first;
}

void GradBandOptimizer::optimize()
{
    cout << "----------------begin optimization-------------" << endl;

    // The first and last points remain unchanged
    // Other intermediate points have three forces: Two from neighbor points and one from gradient
    // F1 = Q_i+1 - Q_i, F2 = Q_i-1 - Q_i, F3 = -(p-p0)*Grad or 0
    // F = F1 + F2 + F3, delta_q = alpha * F
    // Q += delta_q
    // alpha *= beta

    // calculate delta_q for every point, which is a vector<Vector3d> of size (points.rows-2)
    vector<Eigen::Vector3d> delta_qs;
    for (int i = 0; i < points.rows(); ++i)
    {
        // first and last point unmoved
        if (i == 0 || i == points.rows() - 1)
            continue;

        Eigen::Vector3d q1 = points.row(i - 1);
        Eigen::Vector3d q2 = points.row(i);
        Eigen::Vector3d q3 = points.row(i + 1);

        Eigen::Vector3d f1 = q3 - q2;
        Eigen::Vector3d f2 = q1 - q2;

        Eigen::Vector3d grad;
        double dist;
        getDistanceAndGradient(q2, dist, grad);

        Eigen::Vector3d f3;
        if (dist > dist0)
            f3 = Eigen::Vector3d::Zero();
        else
            f3 = -lamda * (dist - dist0) * grad;

        Eigen::Vector3d fi = f1 + f2 + f3;
        Eigen::Vector3d dq = alp * fi;
        delta_qs.push_back(dq);

        // cout << i
        //      << "point \n"
             //      << "q1:\n"
             //      << q1.transpose() << "\nq2:\n"
             //      << q2.transpose() << "\nq3:\n"
             //      << q3.transpose() << "\nf1:\n"
             //      << f1.transpose() << "\nf2:\n"
             //      << f2.transpose() << "\nf3:\n"
             //      << f3.transpose() << "\nfi:\n"
            //  << fi.transpose() << "\ndq:\n"
            //  << dq.transpose() << endl;
    }

    // apply delta q to every points
    for (int i = 0; i < points.rows(); ++i)
    {
        if (i == 0 || i == points.rows() - 1)
            continue;

        points.row(i) += delta_qs[i - 1].transpose();
    }

    // cout << "alp:" << alp << endl;
    alp *= beta;
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