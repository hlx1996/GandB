# Safe trajectory generation in very dense environment based on gradient information and bspline

## Motivation:
Safety is the most important issue when having a uav to do some tasks. 
Corridor based methods generate safe corridors and force the final trajectory to be constrained within them. In this process, a safety margin is chosen according to experience. Large margin could be too conservative and make the optimization problem difficult to solve, while a small one make the resulting trajectory very close to obstacles at some points, which threatens the uav.   
Realizing this drawback, gradient-based methods are proposed which utilize gradient information and push the trajectory far away enough from obstacles. However, these methods also depend heavily on the tuning of parameters. If more weight is put on safety, the resulting trajectory is safe but not smooth, while on the contrary, it is smooth but come too close to obastale. Besides, the trajectory becomes easy to get stuck in local minima as the density of obstacle increase and as a result, the trajectory is not smooth.                       
Thus, we aim to develope a method that take advantage of gradient information and in the mean time do not lose smoothness.

## Problem desciption:
We aim to solve the pose to pose planning problem. During the flying, the uav should always try to fly far enough from the obstacles, at least as far as a pre-defined *margin*. In densely occupied environment, there may be some very narrow passage and there is no trajectory that can pass it without coming into the margin region. At this point, the trajectory should be as far as possible from the obstacles.

## Unclear:
* How to formulate the margin constrain and farest constrain

## Clear:
* bspline basic.
* Minimum snap(jerk) in bspline form.
* ESDF tool and gradient descent of point.
