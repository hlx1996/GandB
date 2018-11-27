# Safe trajectory generation and fast replanning in densely cluttered environment based on gradient information and bspline

## Motivation
Safety is the most important issue when having a UAV to do some tasks. To ensure safety, we need not only reliable global pose to pose planning, but also fast local replanning when some previously undetected obstacles emerge.  

The existing solutions for global pose to pose planning can be roughly categorized into corridor based and gradient based methods. Corridor based methods generate safe corridors and force the final trajectory to be constrained within them. In this process, a safety margin is chosen according to experience. Large margin could result in too conservative trajectory sometimes and make the optimization problem difficult to solve, while a small one make the resulting trajectory very close to obstacles at some points, which threatens the UAV. Realizing this drawback, gradient-based methods are proposed which utilize gradient information and push the trajectory far away enough from obstacles. However, these methods also depend heavily on the tuning of parameters. If more weight is put on safety, the resulting trajectory is safe but not smooth, while on the contrary, it is smooth but come too close to obastale. Besides, the trajectory becomes easy to get stuck in local minima as the density of obstacle increase and as a result, the trajectory is often suboptimal and not smooth either.

In the above mentioned methods, the trajectory is formulated as either piecewise polynomials or piecewise bezier curve. However, I would like to argue that a b-spline representation is more suitable for the purpose of safe trajectory generation and fast replanning:

* Firstly, b-spline is more controlable and predictable. B-spline has the property of strong convex hull. Any b-spline is defined by some control points and can be devided into different segments. The shape of each segment can be controlled by only a small number of related control points. Each segment is always bounded within the convex hull of a few control points. Compared with bezier curve, b-spline get closer to its control polyline(which is formed by connecting the control points with straight line subsequently), which makes its shape more predictable.

* Secondly, b-spline can be modified locally. The change of each control point only affects the shape of the curve locally. This property is desired for fast replanning since we can change a small partition of the trajectory without affecting the global trajectory.

* Thirdly, b-spline ensure fast optimization and derivatives continuity. By the nature of b-spline, a _p_ order b-spline is C^(p-1) continuous. Unlike piecewise polynomials and bezier curve, there is no need to set boundary constrains(derivatives continuity). The only thing to be optimized is the position of each control point.

In conclusion, currently I aim to develope a method that generate safe and smooth trajectory taking advantage of gradient information and b-spline and overcoming the aforementioned drawbacks of previous methods. 

After that I will foucs on the fast replanning problem, especially considering a UAV flying fast in unknown and densely cluttered environment.

## Problem desciption
We aim to solve the pose to pose planning problem. During the flying, the uav should always try to fly far enough from the obstacles, at least as far as a pre-defined *margin*. In densely occupied environment, there may be some very narrow passage and there is no trajectory that can pass it without coming into the margin region. At this point, the trajectory should be as far as possible from the obstacles.

## Problem
* How to formulate the margin constrain and _as far as possible_ constrain.
* What is the ESDF(distance field) of densely cluttered environment look like, how to utilize it.

## Plan
* Spilt and merge of points
 - SQP is very sensitive to large amount of variable (n^3)
* Kinodynamic search


## Result and Plan Log
 - SQP, around from 20-100 ms, result is good, but very slow if number of points is large(50 points, 100ms)
 - constrained slow down SQP (usually double the time)
 - LBFGS is very fast and almost always show good result. Besides it does not get much slower as the number of points increase(50 pts, 1ms; about 100pts, 2ms).


 SDF:
 - sdf_tools: 1000ms for 200x200x50
 - sample function:  for 128x128x128

 ### hybrid a star:

 - when new neighbor localize in the same grid, we need to compare _fscore_ instead of _gscore_
 - no need to constrain neighbor state in neighbor grid (major difference with a star) 
 - currently exclude new state in the same grid


### conversion between hybrid astar trajectory and b-spline 

<!-- - The interval of the b-spline should first selected -> ts
- divide every segment of the b-spline into N parts, each part has the time length of tm = ts/N -> N
- sample at j x tm and get some sample points from the trajectory  
- use this points for least squares problem, solved using Eigen
- when N = 4, i.e., we have 5 sample points for each b-spline segment, the matrix of control point is
[[ 0.00833333333333333       , 0.216666666666667,  0.55, + 0.216666666666667, 0.00833333333333333, 0.0]]
[[ 0.0019775390625           , 0.124910481770833,  0.519645182291667        , 0.328076171875     , + 0.0253824869791667, 0.0]]
[[ 0.000260416666666667      , 0.06171875       ,  0.438020833333333        , 0.438020833333333  , + 0.06171875,  0.000260416666666667]]
[[ 0, + 0.0253824869791667   , 0.328076171875   ,  0.519645182291667        , 0.124910481770833  , + 0.0019775390625]]
[[ 0.0,   0.00833333333333331, 0.216666666666667,  0.55,   0.216666666666667, 0.00833333333333333]]

emmm...this does not work.. because you use a high order polynomial to fit a linear function and a piece-constant function. So the result is ill-conditioning. -->

<!-- - using jerk as input in the path searching stage, and use the same time 
Use jerk as input is not a very good idea. The primitives is easy to drop into the same grid. We can solve this by
enlarging the input jerk and time, but one step will be very long then. This can also be seen from the state transition function(x1 = x0 + v0*t + 0.5*t^2 + 0.16666*j*t^3 and jm*tm <= am), as we need the resulting state to scatter around, either j or t need to be larger. If t is larger, one step is longer.
- Another solution is to lessen the grid size, but this 
slow down the algorithm and need longer time to calculate distance field. -->

- convert the monomial polynomial to b-spline, ok, this is very fast, ok 
  <!-- - but it seems better to use equal distance samples instead of equal time samples, since equal time sample sometimes result in exceeding acceleration(why?). No!!! equal distance sample get very strange result  -->
- consevativeness is a big issue. Maybe we should use lower order bspline so that the control point is closer to real value.
- one shot can be improved

w. burgard


 ### Plan

 - replace signed distance field using the sample function method
  - use ewok's code
 - try 3 order b-spline
 
 - hybrid a star in sdf