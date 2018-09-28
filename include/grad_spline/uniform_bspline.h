#ifndef _UNIFORM_BSPLINE_H_
#define _UNIFORM_BSPLINE_H_

#include <Eigen/Eigen>
#include <iostream>
using namespace std;

class UniformBspline
{
  private:
    /* data */
    int p, n, m;
    Eigen::MatrixXd control_points;
    std::vector<Eigen::MatrixXd> M;  // 012,345
    Eigen::VectorXd u;

    Eigen::VectorXd getU(double u);
    Eigen::MatrixXd getPi(int idx);

  public:
    UniformBspline(Eigen::MatrixXd points, int order, bool auto_extend = true);
    ~UniformBspline();

    void getRegion(double& um, double& um_p);
    Eigen::Vector3d evaluate(double u);
    UniformBspline getDerivative();
};

// control points is a (n+1)x3 matrix
UniformBspline::UniformBspline(Eigen::MatrixXd points, int order, bool auto_extend)
{
    this->p = order;
    if (auto_extend)
    {
        control_points = Eigen::MatrixXd::Zero(points.rows() + 2 * this->p, 3);
        for (int i = 0; i < this->p; ++i)
        {
            control_points.row(i) = points.row(0);
            control_points.row(control_points.rows() - 1 - i) = points.row(points.rows() - 1);
        }
        control_points.block(this->p, 0, points.rows(), 3) = points;
        this->n = points.rows() + 2 * this->p - 1;
    }
    else
    {
        control_points = points;
        this->n = points.rows() - 1;
    }

    this->m = this->n + this->p + 1;

    // calculate knots vector
    this->u = Eigen::VectorXd::Zero(this->m + 1);
    for (int i = 0; i <= this->m; ++i)
    {
        if (i <= this->p)
            this->u(i) = double(-this->p + i);

        else if (i > this->p && i <= this->m - this->p)
        {
            this->u(i) = this->u(i - 1) + 2.0;
        }
        else if (i > this->m - this->p)
        {
            this->u(i) = this->u(i - 1) + 1.0;
        }
    }

    // initialize the M3-5 matrix
    this->M.resize(3);
    Eigen::MatrixXd M3 = Eigen::MatrixXd::Zero(3, 3);
    Eigen::MatrixXd M4 = Eigen::MatrixXd::Zero(4, 4);
    Eigen::MatrixXd M5 = Eigen::MatrixXd::Zero(5, 5);

    M3 << 1.0, 1.0, 0.0, -2.0, 2.0, 0.0, 1.0, -2.0, 1.0;
    M4 << 1.0, 4.0, 1.0, 0.0, -3.0, 0.0, 3.0, 0.0, 3.0, -6.0, 3.0, 0.0, -1.0, 3.0, -3.0, 1.0;
    M5 << 1.0, 11.0, 11.0, 1.0, 0, -4.0, -12.0, 12.0, 4.0, 0, 6.0, -6.0, -6.0, 6.0, 0, -4.0, 12.0, -12.0, 4.0, 0, 1.0,
        -4.0, 6.0, -4.0, 1.0;
    M3 /= 2.0;
    M4 /= 3.0 * 2;
    M5 /= 4.0 * 3.0 * 2.0;
    M[0] = M3;
    M[1] = M4;
    M[2] = M5;

    // show the result
    cout << "p: " << p << "  n: " << n << "  m: " << m << endl;
    cout << "control pts:\n" << control_points << "\nknots:\n" << this->u.transpose() << endl;
    cout << "M3:\n" << M[0] << "\nM4:\n" << M[1] << "\nM5:\n" << M[2] << endl;
}

UniformBspline::~UniformBspline()
{
}

void UniformBspline::getRegion(double& um, double& um_p)
{
    um = this->u(this->p);
    um_p = this->u(this->m - this->p);
}

Eigen::VectorXd UniformBspline::getU(double u)
{
    Eigen::VectorXd uv = Eigen::VectorXd::Zero(this->p + 1);
    uv(0) = 1.0;
    for (int i = 1; i <= this->p; ++i)
    {
        uv(i) = uv(i - 1) * u;
    }
    return uv;
}

Eigen::MatrixXd UniformBspline::getPi(int idx)
{
    Eigen::MatrixXd pi = control_points.block(idx - p, 0, p + 1, 3);
    return pi;
}

Eigen::Vector3d UniformBspline::evaluate(double u)
{
    if (u < this->u(this->p) || u > this->u(this->m - this->p))
        return Eigen::Vector3d::Zero(3);
    // determine which [ui,ui+1] lay in
    int idx = this->p;
    while (true)
    {
        if (this->u(idx + 1) >= u)
            break;
        ++idx;
    }

    u = (u - this->u(idx)) / (this->u(idx + 1) - this->u(idx));

    // get u vector and pi-p -> pi
    Eigen::VectorXd uv = this->getU(u);
    Eigen::MatrixXd pi = this->getPi(idx);
    // cout << "uv:" << uv.transpose() << "\npi: " << pi.transpose() << endl;

    // use p = u'*M*pi
    Eigen::Vector3d val = (uv.transpose() * M[p - 2] * pi).transpose();
    return val;
}

UniformBspline UniformBspline::getDerivative()
{
    // The derivative of a b-spline is also a b-spline, its order become p-1
    // control point Qi = p*(Pi+1-Pi)/(ui+p+1-ui+1)
    Eigen::MatrixXd ctp = Eigen::MatrixXd::Zero(control_points.rows() - 1, 3);
    for (int i = 0; i < ctp.rows(); ++i)
    {
        ctp.row(i) = p * (control_points.row(i + 1) - control_points.row(i)) / (u(i + p + 1) - u(i + 1));
    }
    UniformBspline derivative = UniformBspline(ctp, p - 1, false);
    return derivative;
}

#endif