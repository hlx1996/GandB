#include "grad_spline/grad_band_optimizer.h"
#include <nlopt.hpp>
using namespace std;

class ConstrainData
{
  private:
    /* data */
  public:
    ConstrainData(Eigen::Vector3d fp, Eigen::Vector3d lp, int pn, int vn, int id, int a, int s)
      : first_pt(fp), last_pt(lp), point_num(pn), var_num(vn), idx(id), axis(a), sign(s)
    {
        // show();
    }
    ~ConstrainData(){};

    void show()
    {
        cout << "cons data:\n"
             << first_pt.transpose() << ", " << last_pt.transpose() << "\n"
             << point_num << ", " << var_num << ", " << idx << ", " << axis << ", " << sign << endl;
    }

    Eigen::Vector3d first_pt, last_pt;
    int point_num, var_num;
    int idx, axis, sign;
};

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

void GradBandOptimizer::setParameter(double alpha, double beta, double lamda1, double lamda2, double lamda3, int pow1,
                                     int pow2, double dist0, double scale, int point_opti_num, int algorithm,
                                     double max_vel, double max_acc, double interval)
{
    this->alpha = alpha;
    this->alp = alpha;
    this->beta = beta;
    this->lamda1 = lamda1;
    this->lamda2 = lamda2;
    this->lamda3 = lamda3;
    this->pow1 = pow1;
    this->pow2 = pow2;

    this->dist0 = dist0;

    this->scale = scale;
    this->point_opti_num = point_opti_num;
    this->algorithm = algorithm;

    this->max_vel = max_vel;
    this->max_acc = max_acc;
    this->interval = interval;
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
            f3 = -lamda2 * (dist - dist0) * grad;

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

    // The objective function is f = scale (f1 + lamda2*f2), with f1 = ||q1+q3-2q2||^2
    double f1 = (q1 + q3 - 2 * q2).squaredNorm();

    // f2 = (d-d0)^2, if d<d0; or 0, if d>=d0
    Eigen::Vector3d dist_grad;
    double dist;
    opt->getDistanceAndGradient(q2, dist, dist_grad);

    double f2 = 0.0;
    if (dist < opt->dist0)
        f2 = (dist - opt->dist0) * (dist - opt->dist0);

    double f = opt->scale * (f1 + opt->lamda2 * f2);

    // Then calculate the gradient of the function
    // grad(f1) = -4(q1+q3-2q2)
    Eigen::Vector3d gf1 = -4 * (q1 + q3 - 2 * q2);

    // grad(f2) = 2(d-d0)*grad(dist), if d<d0
    Eigen::Vector3d gf2 = Eigen::Vector3d::Zero();
    if (dist < opt->dist0)
        gf2 = 2 * (dist - opt->dist0) * dist_grad;

    Eigen::Vector3d gf = opt->scale * (gf1 + opt->lamda2 * gf2);

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

// best algorithm is 40: SLSQP(constrained), 11 LBFGS(unconstrained barrier method)
void GradBandOptimizer::optimize3()
{
    this->min_cost = 100000.0;
    nlopt::opt opt(nlopt::algorithm(algorithm), 3 * (points.rows() - 2));
    opt.set_min_objective(GradBandOptimizer::costFunc5, this);

    // add constrains... there are so many, but each constrains is quite simple
    vector<ConstrainData*> cons;
    for (int idx = 0; idx <= this->points.rows() - 2; ++idx)
    {
        for (int axis = 0; axis <= 2; ++axis)
        {
            for (int sign = -1; sign <= 1; ++sign)
            {
                if (sign == 0)
                    continue;

                // this->current_optimize_id = idx;
                // this->current_axis = axis;
                // this->current_sign = sign;
                ConstrainData* con =
                    new ConstrainData(points.row(0).transpose(), points.row(points.rows() - 1).transpose(),
                                      points.rows(), 3 * (points.rows() - 2), idx, axis, sign);
                cons.push_back(con);
                // int con_num = cons.size();

                // opt.add_inequality_constraint(GradBandOptimizer::velConstraint, con, 1e-4);
                // if (idx <= this->points.rows() - 3)
                //     opt.add_inequality_constraint(GradBandOptimizer::accConstraint, con, 1e-4);
            }
        }
    }
    opt.set_maxeval(point_opti_num);
    // opt.set_maxtime(1e-2);
    // opt.set_xtol_rel(1e-4);

    // optimization variables, (x1,y1,z1, x2,y2,z2 ... xn,yn,zn), give initial value
    this->var_num = 3 * (points.rows() - 2);
    vector<double> q(this->var_num);
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
        // optimize and get result
        cout << "begin optimization-------------" << endl;
        nlopt::result result = opt.optimize(q, minf);
        cout << "after all min cost:" << min_cost << " size:" << min_var.size() << endl;
        for (int i = 0; i < points.rows() - 1; ++i)
        {
            if (i == 0 || i == points.rows() - 1)
                continue;

            // points(i, 0) = q[0 + 3 * (i - 1)];
            // points(i, 1) = q[1 + 3 * (i - 1)];
            // points(i, 2) = q[2 + 3 * (i - 1)];
            points(i, 0) = this->min_var[0 + 3 * (i - 1)];
            points(i, 1) = this->min_var[1 + 3 * (i - 1)];
            points(i, 2) = this->min_var[2 + 3 * (i - 1)];
        }

        // cout << "after:  " << points.row(i) << endl;
    }
    catch (std::exception& e)
    {
        // std::cout << "nlopt failed: " << e.what() << std::endl;
    }
}

// this cost function minimize the length; in ensence it optimize the first order derivative
double GradBandOptimizer::costFunc3(const std::vector<double>& x, std::vector<double>& grad, void* func_data)
{
    GradBandOptimizer* opt = reinterpret_cast<GradBandOptimizer*>(func_data);
    grad.resize(opt->var_num);
    for (int i = 0; i < grad.size(); ++i)
    {
        grad[i] = 0;
    }

    // add norm them for the 1->(n-1) terms
    // here we should get the current optimizing point and use its back and front point
    static int optnum = 0;
    ++optnum;
    double f = 0.0, sumf1 = 0.0, sumf2 = 0.0;

    for (int i = 0; i < opt->points.rows(); ++i)
    {
        if (i == 0)
            continue;

        // get the current, back and front points
        Eigen::Vector3d q1, q2;
        if (i == 1)
            q1 = opt->points.row(0);
        else
        {
            q1(0) = x[0 + 3 * (i - 2)];
            q1(1) = x[1 + 3 * (i - 2)];
            q1(2) = x[2 + 3 * (i - 2)];
        }

        if (i == opt->points.rows() - 1)
            q2 = opt->points.row(opt->points.rows() - 1);
        else
        {
            q2(0) = x[0 + 3 * (i - 1)];
            q2(1) = x[1 + 3 * (i - 1)];
            q2(2) = x[2 + 3 * (i - 1)];
        }

        // The objective function is f = scale (f1 + lamda2*f2), with f1 = ||q1-q2||^2
        double f1 = opt->lamda1 * (q1 - q2).squaredNorm();
        sumf1 += f1;

        // f2 = (d-d0)^2, if d<d0; or 0, if d>=d0
        Eigen::Vector3d dist_grad;
        double dist;
        opt->getDistanceAndGradient(q2, dist, dist_grad);

        double f2 = 0.0;
        if (dist < opt->dist0 && i != opt->points.rows() - 1)
            f2 = opt->lamda2 * pow(dist - opt->dist0, opt->pow2);
        sumf2 += f2;

        f += (f1 + f2);

        // Then calculate the gradient of the function
        // for f1 = ||q1+q3-2q2||^2, we much calculate gradient for q1, q2, q3 respectively

        Eigen::Vector3d g1, g2;
        if (i == 1)
            g1 = Eigen::Vector3d::Zero();
        else
            g1 = opt->lamda1 * 2 * (q1 - q2);

        if (i == opt->points.rows() - 1)
            g2 = Eigen::Vector3d::Zero();
        else
            g2 = opt->lamda1 * 2 * (q2 - q1);

        // grad(f2) = 2(d-d0)*grad(dist), if d<d0
        Eigen::Vector3d gf2 = Eigen::Vector3d::Zero();
        if (dist < opt->dist0 && i != opt->points.rows() - 1)
            gf2 = opt->lamda2 * opt->pow2 * pow(dist - opt->dist0, opt->pow2 - 1) * dist_grad;

        // return gradient of function and cost
        if (i != 1)
        {
            grad[0 + 3 * (i - 2)] += g1(0);
            grad[1 + 3 * (i - 2)] += g1(1);
            grad[2 + 3 * (i - 2)] += g1(2);
        }

        if (i != opt->points.rows() - 1)
        {
            grad[0 + 3 * (i - 1)] += g2(0) + gf2(0);
            grad[1 + 3 * (i - 1)] += g2(1) + gf2(1);
            grad[2 + 3 * (i - 1)] += g2(2) + gf2(2);
        }

        // print intermediate result
        // cout << "id: " << opt->current_optimize_id << "\n 3 points:\n"
        //      << q1.transpose() << "\n"
        //      << q2.transpose() << "\n"
        //      << q3.transpose() << "\n f1: " << f1 << "\n f2: " << f2 << " \n gf1: " << gf1.transpose()
        //      << " \n gf2: " << gf2.transpose() << "\n f: " << f << "  grad: " << gf.transpose() << endl;
    }
    cout << optnum << " cost smooth: " << sumf1 << " , grad: " << sumf2 << ", total: " << f << endl;

    // save the min cost result
    if (f < opt->min_cost)
    {
        opt->min_cost = f;
        opt->min_var = x;
    }

    return f;
}

// this cost optimize the curvature; in ensence it optimize the second order derivative
double GradBandOptimizer::costFunc4(const std::vector<double>& x, std::vector<double>& grad, void* func_data)
{
    GradBandOptimizer* opt = reinterpret_cast<GradBandOptimizer*>(func_data);
    grad.resize(opt->var_num);
    for (int i = 0; i < grad.size(); ++i)
    {
        grad[i] = 0;
    }

    // add norm them for the 1->(n-2) terms
    // here we should get the current optimizing point and use its back and front point
    static int optnum = 0;
    ++optnum;
    double f = 0.0, sumf1 = 0.0, sumf2 = 0.0;

    for (int i = 0; i < opt->points.rows(); ++i)
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

        // The objective function is f = scale (f1 + lamda2*f2), with f1 = ||q1+q3-2q2||^2
        double f1 = opt->lamda1 * (q1 + q3 - 2 * q2).squaredNorm();
        sumf1 += f1;

        // f2 = (d-d0)^2, if d<d0; or 0, if d>=d0
        Eigen::Vector3d dist_grad;
        double dist;
        opt->getDistanceAndGradient(q2, dist, dist_grad);

        double f2 = 0.0;
        if (dist < opt->dist0 && i != opt->points.rows() - 1)
            f2 = opt->lamda2 * pow(dist - opt->dist0, opt->pow2);
        sumf2 += f2;

        f += (f1 + f2);

        // Then calculate the gradient of the function
        // for f1 = ||q1+q3-2q2||^2, we much calculate gradient for q1, q2, q3 respectively

        Eigen::Vector3d g1, g2, g3;
        if (i == 1)
            g1 = Eigen::Vector3d::Zero();
        else
            g1 = opt->lamda1 * 2 * (q1 + q3 - 2 * q2);

        g2 = -opt->lamda1 * 4 * (q1 + q3 - 2 * q2);

        if (i == opt->points.rows() - 2)
            g3 = Eigen::Vector3d::Zero();
        else
            g3 = opt->lamda1 * 2 * (q1 + q3 - 2 * q2);

        // grad(f2) = 2(d-d0)*grad(dist), if d<d0
        Eigen::Vector3d gf2 = Eigen::Vector3d::Zero();
        if (dist < opt->dist0 && i != opt->points.rows() - 1)
            gf2 = opt->lamda2 * opt->pow2 * pow(dist - opt->dist0, opt->pow2 - 1) * dist_grad;

        // return gradient of function and cost
        if (i != 1)
        {
            grad[0 + 3 * (i - 2)] += g1(0);
            grad[1 + 3 * (i - 2)] += g1(1);
            grad[2 + 3 * (i - 2)] += g1(2);
        }

        grad[0 + 3 * (i - 1)] += g2(0) + gf2(0);
        grad[1 + 3 * (i - 1)] += g2(1) + gf2(1);
        grad[2 + 3 * (i - 1)] += g2(2) + gf2(2);

        if (i != opt->points.rows() - 2)
        {
            grad[0 + 3 * i] += g3(0);
            grad[1 + 3 * i] += g3(1);
            grad[2 + 3 * i] += g3(2);
        }

        // print intermediate result
        // cout << "id: " << opt->current_optimize_id << "\n 3 points:\n"
        //      << q1.transpose() << "\n"
        //      << q2.transpose() << "\n"
        //      << q3.transpose() << "\n f1: " << f1 << "\n f2: " << f2 << " \n gf1: " << gf1.transpose()
        //      << " \n gf2: " << gf2.transpose() << "\n f: " << f << "  grad: " << gf.transpose() << endl;
    }
    cout << optnum << " cost smooth: " << sumf1 << " , grad: " << sumf2 << ", total: " << f << endl;

    // save the min cost result
    if (f < opt->min_cost)
    {
        opt->min_cost = f;
        opt->min_var = x;
    }

    return f;
}

// this cost optimize the curvature; in ensence it optimize the second order derivative
// it also apply the barrier function method to ensure feasibility
double GradBandOptimizer::costFunc5(const std::vector<double>& x, std::vector<double>& grad, void* func_data)
{
    GradBandOptimizer* opt = reinterpret_cast<GradBandOptimizer*>(func_data);
    grad.resize(opt->var_num);
    for (int i = 0; i < grad.size(); ++i)
    {
        grad[i] = 0;
    }

    // add norm them for the 1->(n-2) terms
    // here we should get the current optimizing point and use its back and front point
    static int optnum = 0;
    ++optnum;
    double f = 0.0, sumf1 = 0.0, sumf2 = 0.0, sumg = 0.0;

    for (int i = 0; i < opt->points.rows(); ++i)
    {
        if (i == 0 || i == opt->points.rows() - 1)
            continue;

        // get the current, back and front points
        Eigen::Vector3d q1, q2, q3;  // 1 previous, 2 present, 3 next point
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

        // The objective function is f = l1*f1 + l2*f2 + r(gv+ga), with f1 = ||q1+q3-2q2||^2
        double f1 = opt->lamda1 * (q1 + q3 - 2 * q2).squaredNorm();
        sumf1 += f1;

        // add the feasibility velocity terms gv, -1/[(mi - mi-1)^2-(vm*dt)^2]
        static double vm = opt->max_vel, am = opt->max_acc, dt = opt->interval;  // limits of vel, acc; and interval
        static double vmdts = vm * dt * vm * dt;                                 // (vm*dt)^2
        Eigen::Vector3d q21 = q2 - q1;                                           // qi - qi-1
        Eigen::VectorXd costv(6);  // denominator of x, y, z; if i = n-2, the last 3 elements are used, else they're
                                   // empty
        double gv = 0.0;

        for (int j = 0; j < 3; ++j)
        {
            double deno = q21(j) * q21(j) - vmdts;
            costv(j) = -1 / deno;
            gv += costv(j);
        }
        // take care if i = n-2
        Eigen::Vector3d q32 = q3 - q2;
        if (i == opt->points.rows() - 2)
        {
            for (int j = 0; j < 3; ++j)
            {
                double deno = q32(j) * q32(j) - vmdts;
                costv(j + 3) = -1 / deno;
                gv += costv(j + 3);
            }
        }
        gv *= opt->lamda3;

        // add the acceleration terms ga, -1/[(mi+1 + mi-1 - 2mi)^2 - (am*dt2)^2]
        static double amdt2s = am * dt * dt * am * dt * dt;  // (am*dt2)^2
        Eigen::Vector3d q132 = q1 + q3 - 2 * q2;
        Eigen::Vector3d costa(3);
        double ga = 0.0;
        for (int j = 0; j < 3; ++j)
        {
            double deno = q132(j) * q132(j) - amdt2s;
            costa(j) = -1 / deno;
            ga += costa(j);
        }
        ga *= opt->lamda3;

        // show the points relation and g cost
        for (int j = 0; j < 3; ++j)
        {
            // if (costv(j) < 0 || costa(j) < 0)
            //     cout << "q21:" << q21.transpose() << "\nv cost:" << costv.head(3).transpose()
            //          << "\nq132:" << q132.transpose() << "\na cost:" << costa.transpose() << "\n\n";
        }

        sumg += ga + gv;

        // f2 = (d-d0)^2, if d<d0; or 0, if d>=d0
        Eigen::Vector3d dist_grad;
        double dist;
        opt->getDistanceAndGradient(q2, dist, dist_grad);

        double f2 = 0.0;
        if (dist < opt->dist0 && i != opt->points.rows() - 1)
            f2 = opt->lamda2 * pow(dist - opt->dist0, opt->pow2);
        sumf2 += f2;

        f += (f1 + f2 + gv + ga);

        // Then calculate the gradient of the function
        // for f1 = ||q1+q3-2q2||^2, we much calculate gradient for q1, q2, q3 respectively
        Eigen::Vector3d g1, g2, g3;
        if (i == 1)
            g1 = Eigen::Vector3d::Zero();
        else
            g1 = opt->lamda1 * 2 * q132;

        g2 = -opt->lamda1 * 4 * q132;

        if (i == opt->points.rows() - 2)
            g3 = Eigen::Vector3d::Zero();
        else
            g3 = opt->lamda1 * 2 * q132;

        // add the gradient of gv and ga
        for (int j = 0; j < 3; ++j)
        {
            // for gv, grad(B/mi) = 1/[]^2 * 2(mi - mi-1)
            grad[j + 3 * (i - 1)] += opt->lamda3 * costv(j) * costv(j) * 2 * q21(j);

            // for gv, grad(B/mi-1) = 1/[]^2 * 2(mi-1 - mi)
            if (i != 1)
                grad[j + 3 * (i - 2)] += -opt->lamda3 * costv(j) * costv(j) * 2 * q21(j);

            // take care when i = n-2
            if (i == opt->points.rows() - 2)
            {
                grad[j + 3 * (i - 1)] += -opt->lamda3 * costv(j + 3) * costv(j + 3) * 2 * q32(j);
            }

            // then for ga, grad(B/mi) = 1/[]^2 *(-4)*(mi+1 + mi-1 -2mi2)
            grad[j + 3 * (i - 1)] += opt->lamda3 * costa(j) * costa(j) * (-4) * q132(j);
            if (i != 1)
                grad[j + 3 * (i - 2)] += opt->lamda3 * costa(j) * costa(j) * 2 * q132(j);
            if (i != opt->points.rows() - 2)
                grad[j + 3 * i] += opt->lamda3 * costa(j) * costa(j) * 2 * q132(j);
        }

        // grad(f2) = 2(d-d0)*grad(dist), if d<d0
        Eigen::Vector3d gf2 = Eigen::Vector3d::Zero();
        if (dist < opt->dist0 && i != opt->points.rows() - 1)
            gf2 = opt->lamda2 * opt->pow2 * pow(dist - opt->dist0, opt->pow2 - 1) * dist_grad;

        // return gradient of function and cost
        if (i != 1)
        {
            grad[0 + 3 * (i - 2)] += g1(0);
            grad[1 + 3 * (i - 2)] += g1(1);
            grad[2 + 3 * (i - 2)] += g1(2);
        }

        grad[0 + 3 * (i - 1)] += g2(0) + gf2(0);
        grad[1 + 3 * (i - 1)] += g2(1) + gf2(1);
        grad[2 + 3 * (i - 1)] += g2(2) + gf2(2);

        if (i != opt->points.rows() - 2)
        {
            grad[0 + 3 * i] += g3(0);
            grad[1 + 3 * i] += g3(1);
            grad[2 + 3 * i] += g3(2);
        }

        // print intermediate result
        // cout << "id: " << opt->current_optimize_id << "\n 3 points:\n"
        //      << q1.transpose() << "\n"
        //      << q2.transpose() << "\n"
        //      << q3.transpose() << "\n f1: " << f1 << "\n f2: " << f2 << " \n gf1: " << gf1.transpose()
        //      << " \n gf2: " << gf2.transpose() << "\n f: " << f << "  grad: " << gf.transpose() << endl;
    }
    cout << optnum << " smooth: " << sumf1 << " , grad: " << sumf2 << " , barrier: " << sumg << ", total: " << f
         << "\n\n";

    // save the min cost result
    if (f < opt->min_cost)
    {
        opt->min_cost = f;
        opt->min_var = x;
    }

    return f;
}

double GradBandOptimizer::velConstraint(const std::vector<double>& x, std::vector<double>& grad, void* data)
{
    ConstrainData* cons = reinterpret_cast<ConstrainData*>(data);
    grad.resize(cons->var_num);

    // get the current idx of points
    // int idx = cons->opt->current_optimize_id;
    // int axis = cons->opt->current_axis;
    // int sign = cons->opt->current_sign;
    int idx = cons->idx;
    int axis = cons->axis;
    int sign = cons->sign;
    // cout << "vel cons:" << endl;
    // cons->show();

    // the control points of velocity is (pi+1-pi)/dt, when the interal is dt
    // we try interval = 2
    // optimization variables, (x1,y1,z1, x2,y2,z2 ... xn,yn,zn)
    // notice that the first and last point is not added in the optimized variables
    double mi, mi1, vi;
    for (int i = 0; i < cons->var_num - 1; ++i)
        grad[i] = 0.0;

    // constrain
    if (idx == 0)
    {
        mi = cons->first_pt(axis);
        mi1 = x[axis];

        grad[axis] = double(sign);
    }
    else if (idx == cons->point_num - 2)
    {
        mi = x[3 * (idx - 1) + axis];
        mi1 = cons->last_pt(axis);

        grad[3 * (idx - 1) + axis] = -double(sign);
    }
    else
    {
        mi = x[3 * (idx - 1) + axis];
        mi1 = x[3 * idx + axis];

        grad[3 * (idx - 1) + axis] = -double(sign);
        grad[3 * idx + axis] = double(sign);
    }

    // NLopt use g(x) <= 0 constraint
    static double vm = 2.5, dt = 2.0;
    return sign * (mi1 - mi) - vm * dt;
}

double GradBandOptimizer::accConstraint(const std::vector<double>& x, std::vector<double>& grad, void* data)
{
    ConstrainData* cons = reinterpret_cast<ConstrainData*>(data);
    grad.resize(cons->var_num);

    // get the current idx of points
    int idx = cons->idx;
    int axis = cons->axis;
    int sign = cons->sign;
    // cout << "acc cons" << endl;
    // cons->show();

    // the control points of acceleration is (pi+2 - 2pi+1 + pi)/dt^2, when the interal is dt
    // we try interval = 2
    // optimization variables, (x1,y1,z1, x2,y2,z2 ... xn,yn,zn)
    // notice that the first and last point is not added in the optimized variables
    double mi, mi1, mi2, ai;
    for (int i = 0; i < cons->var_num - 1; ++i)
        grad[i] = 0.0;

    // constrain
    if (idx == 0)
    {
        mi = cons->first_pt(axis);
        mi1 = x[axis];
        mi2 = x[axis + 3];

        grad[axis] = -2.0 * double(sign);
        grad[axis + 3] = double(sign);
    }
    else if (idx == cons->point_num - 3)
    {
        mi = x[3 * (idx - 1) + axis];
        mi1 = x[3 * idx + axis];
        mi2 = cons->last_pt(axis);

        grad[3 * (idx - 1) + axis] = double(sign);
        grad[3 * idx + axis] = -2.0 * double(sign);
    }
    else
    {
        mi = x[3 * (idx - 1) + axis];
        mi1 = x[3 * idx + axis];
        mi2 = x[3 * (idx + 1) + axis];

        grad[3 * (idx - 1) + axis] = double(sign);
        grad[3 * idx + axis] = -2.0 * double(sign);
        grad[3 * (idx + 1) + axis] = double(sign);
    }

    // NLopt use g(x) <= 0 constraint
    static double am = 2.5, dt = 2.0;
    return sign * (mi2 - 2 * mi1 + mi) - am * dt * dt;
}