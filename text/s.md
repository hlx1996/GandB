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
* bspline basic. write code of b-spline
* Minimum snap(jerk) in bspline form.
* ESDF tool and gradient method for moving a point.
* move control point in the gradient of jerk.
