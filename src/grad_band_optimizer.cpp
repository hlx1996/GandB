#include "grad_spline/grad_band_optimizer.h"
#include <nlopt.hpp>
using namespace std;

GradBandOptimizer::GradBandOptimizer(Eigen::MatrixXd points, sdf_tools::SignedDistanceField* sdf, double res)
{
    this->sdf = sdf;
    this->resolution = res;
    this->points = points;

    // cout << "points:\n" << this->points << endl;
}

GradBandOptimizer::~GradBandOptimizer()
{
}

Eigen::MatrixXd GradBandOptimizer::getPoints()
{
    return this->points;
}

void GradBandOptimizer::setParameter(double alpha, double beta, double lamda, double dist0, double scale,
                                     int point_opti_num, int algorithm)
{
    this->alpha = alpha;
    this->alp = alpha;
    this->beta = beta;
    this->lamda = lamda;
    this->dist0 = dist0;

    this->scale = scale;
    this->point_opti_num = point_opti_num;
    this->algorithm = algorithm;
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

        current_optimize_id = i;

        Eigen::Vector3d q1 = points.row(i - 1);
        Eigen::Vector3d q2 = points.row(i);
        Eigen::Vector3d q3 = points.row(i + 1);

        Eigen::Vector3d f1 = q3 - q2;
        Eigen::Vector3d f2 = q1 - q2;

        Eigen::Vector3d grad;
        double dist;
        getDistanceAndGradient(q2, dist, grad);

        Eigen::Vector3d f3;
        if (dist >= dist0)
            f3 = Eigen::Vector3d::Zero();
        else
            f3 = -lamda * (dist - dist0) * grad;

        Eigen::Vector3d fi = f1 + f2 + f3;
        Eigen::Vector3d dq = alp * fi;
        delta_qs.push_back(dq);

        // test applying the dq in sequence ----(1)
        points.row(i) += dq.transpose();

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

    // apply delta q to every points in the same time, on the contrary to (1)
    // for (int i = 0; i < points.rows(); ++i)
    // {
    //     if (i == 0 || i == points.rows() - 1)
    //         continue;

    //     points.row(i) += delta_qs[i - 1].transpose();
    // }

    // cout << "alp:" << alp << endl;
    alp *= beta;
}

// best algorithm is PRECOND-NEWTON
void GradBandOptimizer::optimize2()
{
    nlopt::opt opt(nlopt::algorithm(algorithm), 3);
    opt.set_min_objective(GradBandOptimizer::costFunc2, this);

    opt.set_maxeval(point_opti_num);
    // opt.set_maxtime(1e-2);
    // opt.set_xtol_rel(1e-4);

    for (int i = 0; i < points.rows(); ++i)
    {
        if (i == 0 || i == points.rows() - 1)
            continue;

        current_optimize_id = i;
        // cout << "optimizing the " << i << " th point------------------------------- \n";

        // Here we optimize each point using NLopt solver
        // For each point, create a NLopt solver, set costfunc and constrains
        // Then set stopping criterion and give a initial guess
        // Finally call optimize() with vector<double> as input and output

        vector<double> q(3);
        q[0] = points.row(i)(0);
        q[1] = points.row(i)(1);
        q[2] = points.row(i)(2);
        double minf;

        // cout << "before:  " << points.row(i) << endl;

        try
        {
            nlopt::result result = opt.optimize(q, minf);
            // std::cout << "found minimum at f(" << x[0] << "," << x[1] << ") = " << std::setprecision(10) << minf
            //           << std::endl;
            points.row(i)(0) = q[0];
            points.row(i)(1) = q[1];
            points.row(i)(2) = q[2];

            // cout << "after:  " << points.row(i) << endl;
        }
        catch (std::exception& e)
        {
            // std::cout << "nlopt failed: " << e.what() << std::endl;
        }
    }
}

// best algorithm is 11: LBFGS
void GradBandOptimizer::optimize3()
{
    nlopt::opt opt(nlopt::algorithm(algorithm), 3 * (points.rows() - 2));
    opt.set_min_objective(GradBandOptimizer::costFunc3, this);

    opt.set_maxeval(point_opti_num);
    // opt.set_maxtime(1e-2);
    // opt.set_xtol_rel(1e-4);

    // optimization variables, (x1,y1,z1, x2,y2,z2 ... xn,yn,zn)
    vector<double> q(3 * (points.rows() - 2));
    double minf;
    for (int i = 0; i < points.rows() - 1; ++i)
    {
        if (i == 0 || i == points.rows() - 1)
            continue;

        q[0 + 3 * (i - 1)] = points(i, 0);
        q[1 + 3 * (i - 1)] = points(i, 1);
        q[2 + 3 * (i - 1)] = points(i, 2);
    }

    try
    {
        nlopt::result result = opt.optimize(q, minf);
        for (int i = 0; i < points.rows() - 1; ++i)
        {
            if (i == 0 || i == points.rows() - 1)
                continue;

            points(i, 0) = q[0 + 3 * (i - 1)];
            points(i, 1) = q[1 + 3 * (i - 1)];
            points(i, 2) = q[2 + 3 * (i - 1)];
        }

        // cout << "after:  " << points.row(i) << endl;
    }
    catch (std::exception& e)
    {
        // std::cout << "nlopt failed: " << e.what() << std::endl;
    }
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

double GradBandOptimizer::costFunc2(const std::vector<double>& x, std::vector<double>& grad, void* func_data)
{
    GradBandOptimizer* opt = reinterpret_cast<GradBandOptimizer*>(func_data);

    // here we should get the current optimizing point and use its back and front point
    Eigen::Vector3d q1 = opt->points.row(opt->current_optimize_id - 1);
    Eigen::Vector3d q3 = opt->points.row(opt->current_optimize_id + 1);
    Eigen::Vector3d q2;
    q2(0) = x[0];
    q2(1) = x[1];
    q2(2) = x[2];

    // The objective function is f = scale (f1 + lamda*f2), with f1 = ||q1+q3-2q2||^2
    double f1 = (q1 + q3 - 2 * q2).squaredNorm();

    // f2 = (d-d0)^2, if d<d0; or 0, if d>=d0
    Eigen::Vector3d dist_grad;
    double dist;
    opt->getDistanceAndGradient(q2, dist, dist_grad);

    double f2 = 0.0;
    if (dist < opt->dist0)
        f2 = (dist - opt->dist0) * (dist - opt->dist0);

    double f = opt->scale * (f1 + opt->lamda * f2);

    // Then calculate the gradient of the function
    // grad(f1) = -4(q1+q3-2q2)
    Eigen::Vector3d gf1 = -4 * (q1 + q3 - 2 * q2);

    // grad(f2) = 2(d-d0)*grad(dist), if d<d0
    Eigen::Vector3d gf2 = Eigen::Vector3d::Zero();
    if (dist < opt->dist0)
        gf2 = 2 * (dist - opt->dist0) * dist_grad;

    Eigen::Vector3d gf = opt->scale * (gf1 + opt->lamda * gf2);

    // return gradient of function and cost
    grad.resize(3);
    grad[0] = gf(0);
    grad[1] = gf(1);
    grad[2] = gf(2);

    // print intermediate result
    // cout << "id: " << opt->current_optimize_id << "\n 3 points:\n"
    //      << q1.transpose() << "\n"
    //      << q2.transpose() << "\n"
    //      << q3.transpose() << "\n f1: " << f1 << "\n f2: " << f2 << " \n gf1: " << gf1.transpose()
    //      << " \n gf2: " << gf2.transpose() << "\n f: " << f << "  grad: " << gf.transpose() << endl;

    return f;
}

double GradBandOptimizer::costFunc3(const std::vector<double>& x, std::vector<double>& grad, void* func_data)
{
    GradBandOptimizer* opt = reinterpret_cast<GradBandOptimizer*>(func_data);
    grad.resize(3 * (opt->points.rows() - 2));

    // here we should get the current optimizing point and use its back and front point
    double f = 0.0;
    for (int i = 0; i < opt->points.rows() - 1; ++i)
    {
        if (i == 0 || i == opt->points.rows() - 1)
            continue;

        // get the current, back and front points
        Eigen::Vector3d q1, q2, q3;
        if (i == 1)
            q1 = opt->points.row(0);
        else
        {
            q1(0) = x[0 + 3 * (i - 2)];
            q1(1) = x[1 + 3 * (i - 2)];
            q1(2) = x[2 + 3 * (i - 2)];
        }

        q2(0) = x[0 + 3 * (i - 1)];
        q2(1) = x[1 + 3 * (i - 1)];
        q2(2) = x[2 + 3 * (i - 1)];

        if (i == opt->points.rows() - 2)
            q3 = opt->points.row(opt->points.rows() - 1);
        else
        {
            q3(0) = x[0 + 3 * i];
            q3(1) = x[1 + 3 * i];
            q3(2) = x[2 + 3 * i];
        }

        // The objective function is f = scale (f1 + lamda*f2), with f1 = ||q1+q3-2q2||^2
        double f1 = (q1 + q3 - 2 * q2).squaredNorm();

        // f2 = (d-d0)^2, if d<d0; or 0, if d>=d0
        Eigen::Vector3d dist_grad;
        double dist;
        opt->getDistanceAndGradient(q2, dist, dist_grad);

        double f2 = 0.0;
        if (dist < opt->dist0)
            f2 = (dist - opt->dist0) * (dist - opt->dist0);

        f += opt->scale * (f1 + opt->lamda * f2);

        // Then calculate the gradient of the function
        // grad(f1) = -4(q1+q3-2q2)
        Eigen::Vector3d gf1 = -4 * (q1 + q3 - 2 * q2);

        // grad(f2) = 2(d-d0)*grad(dist), if d<d0
        Eigen::Vector3d gf2 = Eigen::Vector3d::Zero();
        if (dist < opt->dist0)
            gf2 = 2 * (dist - opt->dist0) * dist_grad;

        Eigen::Vector3d gf = opt->scale * (gf1 + opt->lamda * gf2);

        // return gradient of function and cost
        grad[0 + 3 * (i - 1)] = gf(0);
        grad[1 + 3 * (i - 1)] = gf(1);
        grad[2 + 3 * (i - 1)] = gf(2);

        // print intermediate result
        // cout << "id: " << opt->current_optimize_id << "\n 3 points:\n"
        //      << q1.transpose() << "\n"
        //      << q2.transpose() << "\n"
        //      << q3.transpose() << "\n f1: " << f1 << "\n f2: " << f2 << " \n gf1: " << gf1.transpose()
        //      << " \n gf2: " << gf2.transpose() << "\n f: " << f << "  grad: " << gf.transpose() << endl;
    }
    return f;
}