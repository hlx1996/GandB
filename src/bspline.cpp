#include "grad_spline/uniform_bspline.h"

int main(int argc, char const* argv[])
{
    /* code */
    Eigen::MatrixXd points = Eigen::MatrixXd::Zero(9, 3);
    points << 1.0, 2.0, 2.0, 2.0, 3.0, 2.0, 3.0, 5.0, 2.0, 4.0, 7.0, 2.0, 5.0, 3.0, 2.0, 6.0, 6.0, 2.0, 7.0, 8.0, 2.0,
        8.0, 9.0, 2.0, 9.0, 1.0, 2.0;

    cout << "points:\n" << points << endl;

    UniformBspline bspline = UniformBspline(points, 4, 2.0);
    double a, b;
    bspline.getRegion(a, b);
    cout << "region: " << a << "," << b << endl;

    cout << "evaluate:" << endl;
    cout << bspline.evaluate(3.0) << endl;
    cout << bspline.evaluate(5.0) << endl;
    cout << bspline.evaluate(7.0) << endl;

    UniformBspline d1 = bspline.getDerivative();
    cout << "derivative:" << endl;
    cout << d1.evaluate(3.0) << endl;
    cout << d1.evaluate(5.0) << endl;
    cout << d1.evaluate(7.0) << endl;
    cout << d1.evaluate(20.0) << endl;
    return 0;
}
