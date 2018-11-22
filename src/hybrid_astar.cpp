#include "grad_spline/hybrid_astar.h"
#include <sstream>

using namespace std;
using namespace Eigen;
using namespace sdf_tools;

void gridPathFinder::initGridNodeMap(double _resolution, Vector3d global_xyz_l)
{
    gl_xl = global_xyz_l(0);
    gl_yl = global_xyz_l(1);
    gl_zl = global_xyz_l(2);

    resolution = _resolution;
    inv_resolution = 1.0 / _resolution;

    w_time = 30;

    GridNodeMap = new GridNodePtr**[GLX_SIZE];
    for (int i = 0; i < GLX_SIZE; i++)
    {
        GridNodeMap[i] = new GridNodePtr*[GLY_SIZE];
        for (int j = 0; j < GLY_SIZE; j++)
        {
            GridNodeMap[i][j] = new GridNodePtr[GLZ_SIZE];
            for (int k = 0; k < GLZ_SIZE; k++)
            {
                Vector3i tmpIdx(i, j, k);
                Vector3d pos = gridIndex2coord(tmpIdx);
                GridNodeMap[i][j][k] = new GridNode(tmpIdx, pos);
            }
        }
    }
}

void gridPathFinder::linkLocalMap(CollisionMapGrid* local_map, Vector3d xyz_l)
{
    Vector3d coord;
    for (int64_t i = 0; i < X_SIZE; i++)
    {
        for (int64_t j = 0; j < Y_SIZE; j++)
        {
            for (int64_t k = 0; k < Z_SIZE; k++)
            {
                coord(0) = xyz_l(0) + (double)(i + 0.5) * resolution;
                coord(1) = xyz_l(1) + (double)(j + 0.5) * resolution;
                coord(2) = xyz_l(2) + (double)(k + 0.5) * resolution;

                Vector3i index = coord2gridIndex(coord);

                if (index(0) >= GLX_SIZE || index(1) >= GLY_SIZE || index(2) >= GLZ_SIZE || index(0) < 0 ||
                    index(1) < 0 || index(2) < 0)
                    continue;

                GridNodePtr ptr = GridNodeMap[index(0)][index(1)][index(2)];
                ptr->id = 0;
                ptr->occupancy = local_map->Get(i, j, k).first.occupancy;
            }
        }
    }
}

void gridPathFinder::resetLocalMap()
{
    // ROS_WARN("expandedNodes size : %d", expandedNodes.size());
    for (auto tmpPtr : expandedNodes)
    {
        tmpPtr->occupancy = 0;  // forget the occupancy
        tmpPtr->id = 0;
        tmpPtr->cameFrom = NULL;
        tmpPtr->gScore = inf;
        tmpPtr->fScore = inf;
    }

    for (auto ptr : openSet)
    {
        GridNodePtr tmpPtr = ptr.second;
        tmpPtr->occupancy = 0;  // forget the occupancy
        tmpPtr->id = 0;
        tmpPtr->cameFrom = NULL;
        tmpPtr->gScore = inf;
        tmpPtr->fScore = inf;
    }

    expandedNodes.clear();
    // ROS_WARN("local map reset finish");
}

GridNodePtr gridPathFinder::pos2gridNodePtr(Vector3d pos)
{
    Vector3i idx = coord2gridIndex(pos);
    GridNodePtr grid_ptr = new GridNode(idx, pos);

    Eigen::VectorXd state(6);
    state.head(3) = pos;
    state.tail(3) = Eigen::Vector3d::Zero();
    grid_ptr->state = state;

    return grid_ptr;
}

Vector3d gridPathFinder::gridIndex2coord(Vector3i index)
{
    Vector3d pt;
    // cell_x_size_ * ((double)x_index + 0.5), cell_y_size_ * ((double)y_index +
    // 0.5), cell_z_size_ * ((double)z_index + 0.5)

    pt(0) = ((double)index(0) + 0.5) * resolution + gl_xl;
    pt(1) = ((double)index(1) + 0.5) * resolution + gl_yl;
    pt(2) = ((double)index(2) + 0.5) * resolution + gl_zl;

    /*pt(0) = (double)index(0) * resolution + gl_xl + 0.5 * resolution;
    pt(1) = (double)index(1) * resolution + gl_yl + 0.5 * resolution;
    pt(2) = (double)index(2) * resolution + gl_zl + 0.5 * resolution;*/
    return pt;
}

Vector3i gridPathFinder::coord2gridIndex(Vector3d pt)
{
    Vector3i idx;
    idx << min(max(int((pt(0) - gl_xl) * inv_resolution), 0), GLX_SIZE - 1),
        min(max(int((pt(1) - gl_yl) * inv_resolution), 0), GLY_SIZE - 1),
        min(max(int((pt(2) - gl_zl) * inv_resolution), 0), GLZ_SIZE - 1);

    return idx;
}

double gridPathFinder::getDiagHeu(GridNodePtr node1, GridNodePtr node2)
{
    double dx = abs(node1->index(0) - node2->index(0));
    double dy = abs(node1->index(1) - node2->index(1));
    double dz = abs(node1->index(2) - node2->index(2));

    double h;
    int diag = min(min(dx, dy), dz);
    dx -= diag;
    dy -= diag;
    dz -= diag;

    if (dx == 0)
    {
        h = 1.0 * sqrt(3.0) * diag + sqrt(2.0) * min(dy, dz) + 1.0 * abs(dy - dz);
    }
    if (dy == 0)
    {
        h = 1.0 * sqrt(3.0) * diag + sqrt(2.0) * min(dx, dz) + 1.0 * abs(dx - dz);
    }
    if (dz == 0)
    {
        h = 1.0 * sqrt(3.0) * diag + sqrt(2.0) * min(dx, dy) + 1.0 * abs(dx - dy);
    }
    return h;
}

double gridPathFinder::getManhHeu(GridNodePtr node1, GridNodePtr node2)
{
    double dx = abs(node1->index(0) - node2->index(0));
    double dy = abs(node1->index(1) - node2->index(1));
    double dz = abs(node1->index(2) - node2->index(2));

    return dx + dy + dz;
}

double gridPathFinder::getEuclHeu(GridNodePtr node1, GridNodePtr node2)
{
    return (node2->index - node1->index).norm();
}

double gridPathFinder::getHeu(GridNodePtr node1, GridNodePtr node2)
{
    return tie_breaker * getDiagHeu(node1, node2);
    // return tie_breaker * getEuclHeu(node1, node2);
}

vector<GridNodePtr> gridPathFinder::retrievePath(GridNodePtr current)
{
    vector<GridNodePtr> path;
    path.push_back(current);

    while (current->cameFrom != NULL)
    {
        current = current->cameFrom;
        path.push_back(current);
    }

    return path;
}

vector<GridNodePtr> gridPathFinder::getVisitedNodes()
{
    vector<GridNodePtr> d;
    for (int i = 0; i < GLX_SIZE; i++)
        for (int j = 0; j < GLY_SIZE; j++)
            for (int k = 0; k < GLZ_SIZE; k++)
            {
                if (GridNodeMap[i][j][k]->id != 0)
                    // if(GridNodeMap[i][j][k]->id == -1)
                    d.push_back(GridNodeMap[i][j][k]);
            }

    ROS_WARN("d size : %d", d.size());
    return d;
}

/*bool gridPathFinder::minClearance()
{
    neighborPtr->occupancy > 0.5
}
*/
void gridPathFinder::AstarSearch(Vector3d start_pt, Vector3d start_vel, Vector3d end_pt, Vector3d end_vel)
{
    cout << "[Hybrid A Star]: begin" << endl;
    ros::Time time_1 = ros::Time::now();
    double vis_time = 0.0;

    Vector3i start_idx = coord2gridIndex(start_pt);
    Vector3i end_idx = coord2gridIndex(end_pt);

    GridNodePtr startPtr = GridNodeMap[start_idx(0)][start_idx(1)][start_idx(2)];
    GridNodePtr endPtr = GridNodeMap[end_idx(0)][end_idx(1)][end_idx(2)];

    startPtr->state.head(3) = start_pt;
    startPtr->state.segment(3, 3) = start_vel;
    // startPtr->state.tail(3) = Eigen::Vector3d::Zero(3);
    endPtr->state.head(3) = end_pt;
    endPtr->state.segment(3, 3) = end_vel;
    // endPtr->state.tail(3) = Eigen::Vector3d::Zero(3);

    startPtr->gScore = 0;
    startPtr->fScore = getKinoDynamicHeu(startPtr, endPtr);
    startPtr->id = 1;  // put start node in open set
    startPtr->coord = start_pt;

    openSet.clear();
    openSet.insert(make_pair(startPtr->fScore, startPtr));  // put start in open set

    GridNodePtr neighbor_ptr = NULL;
    GridNodePtr current_node = NULL;
    terminate_ptr = NULL;

    double tentative_gScore;
    int num_iter = 0;

    // num_ope = 0;
    // time_in_forward = 0.0;
    is_shot_succ = false;
    N_max_shot = 10;
    cnt_shot = 0;
    dis_shot = (startPtr->state.head(3) - endPtr->state.head(3)).norm();

    while (!openSet.empty())
    {
        ++num_iter;
        current_node = openSet.begin()->second;

        int difference = 5;
        if (abs(current_node->index(0) - endPtr->index(0)) <= difference &&
            abs(current_node->index(1) - endPtr->index(1)) <= difference &&
            abs(current_node->index(2) - endPtr->index(2)) <= difference)
        {
            // the final distance is reached using one shot
            // shotHeu(current_node, endPtr)
            ROS_WARN("[hybrid Astar]Reach goal..");
            // cout << "goal coord: " << endl << current_node->real_coord << endl;
            cout << "total number of iteration used in hybrid Astar: " << num_iter << endl;
            ros::Time time_2 = ros::Time::now();
            ROS_WARN("Time consume in hybrid A star path finding is %f", (time_2 - time_1).toSec() - vis_time);
            gridPath = retrievePath(current_node);
            terminate_ptr = current_node;
            has_path = true;
            return;
        }
        openSet.erase(openSet.begin());
        current_node->id = -1;  // move node from open to closed set.
        expandedNodes.push_back(current_node);

        // get the state of current node
        Eigen::VectorXd state = current_node->state;
        // cout << "\ncurernt state:" << state.transpose() << endl;

        // get neighbor of this node
        Neighbors neighbors;
        getNeighbor(current_node, endPtr, neighbors, vis_time);
        vector<pair<Vector3i, KinoState>> neighbor_data = neighbors.getData();

        // iterate the neighbors
        for (int i = 0; i < int(neighbor_data.size()); ++i)
        {
            // get the neighbor node
            Eigen::Vector3i diff, neighbor_idx;
            diff = neighbor_data[i].first;
            neighbor_idx(0) = current_node->index(0) + diff(0);
            neighbor_idx(1) = current_node->index(1) + diff(1);
            neighbor_idx(2) = current_node->index(2) + diff(2);

            KinoState neighbor = neighbor_data[i].second;
            neighbor_ptr = GridNodeMap[neighbor_idx(0)][neighbor_idx(1)][neighbor_idx(2)];

            double edge_cost = neighbor.edge_cost;
            double heu = neighbor.heu;
            double optimal_time = neighbor.optimal_time;
            double duration = neighbor.duration;
            Eigen::Vector3d input = neighbor.input;
            Eigen::VectorXd x1 = neighbor.state;

            tentative_gScore = current_node->gScore + edge_cost;
            if (neighbor_ptr->index == current_node->index)  // in the same grid, need compare
            {
                double current_fscore = current_node->fScore;
                double tentative_fscore = current_node->gScore + edge_cost + heu;

                if (tentative_fscore < current_fscore + tie_breaker)
                {
                    // neighborPtr -> input = u;
                    neighbor_ptr->state = x1;
                    neighbor_ptr->gScore = tentative_gScore;
                    neighbor_ptr->fScore = tentative_fscore;
                    neighbor_ptr->optimal_time = optimal_time;
                    neighbor_ptr->duration = duration;
                    neighbor_ptr->input = input;

                    if (neighbor_ptr->id == 1)
                    {
                        if (num_iter < 2) ROS_WARN(" re-insert, take place");

                        openSet.erase(neighbor_ptr->nodeMapIt);
                    }

                    neighbor_ptr->id = 1;
                    neighbor_ptr->nodeMapIt =
                        openSet.insert(make_pair(neighbor_ptr->fScore,
                                                 neighbor_ptr));  // put neighbor in open set and record it.
                }
            }
            else if (neighbor_ptr->id != -1)  // in close set
            {
                if (neighbor_ptr->id != 1)  // discover a new node
                {
                    neighbor_ptr->id = 1;
                    neighbor_ptr->state = x1;
                    neighbor_ptr->optimal_time = optimal_time;
                    neighbor_ptr->duration = duration;
                    neighbor_ptr->input = input;
                    neighbor_ptr->cameFrom = current_node;
                    neighbor_ptr->gScore = tentative_gScore;
                    neighbor_ptr->fScore = neighbor_ptr->gScore + heu;
                    // neighbor_ptr->fScore = neighbor_ptr->gScore + 0.0;
                    neighbor_ptr->nodeMapIt =
                        openSet.insert(make_pair(neighbor_ptr->fScore,
                                                 neighbor_ptr));  // put neighbor in open set and record it.
                    continue;
                }
                else if (tentative_gScore <= neighbor_ptr->gScore)  // already in open set, need compare
                                                                    // and update
                {
                    neighbor_ptr->state = x1;
                    neighbor_ptr->optimal_time = optimal_time;
                    neighbor_ptr->duration = duration;
                    neighbor_ptr->input = input;
                    neighbor_ptr->cameFrom = current_node;
                    neighbor_ptr->gScore = tentative_gScore;
                    // neighbor_ptr->fScore = neighbor_ptr->gScore + 0.0;
                    neighbor_ptr->fScore = tentative_gScore + heu;
                    openSet.erase(neighbor_ptr->nodeMapIt);
                    neighbor_ptr->nodeMapIt =
                        openSet.insert(make_pair(neighbor_ptr->fScore,
                                                 neighbor_ptr));  // put neighbor in open set and record it.
                }
            }
        }
    }

    ros::Time time_2 = ros::Time::now();
    ROS_WARN("End. Time consume in hybrid A star path finding is %f", (time_2 - time_1).toSec());
}

vector<Vector3d> gridPathFinder::getPath()
{
    vector<Vector3d> path;

    for (auto ptr : gridPath) path.push_back(ptr->coord);

    reverse(path.begin(), path.end());
    return path;
}

std::vector<GridNodePtr> gridPathFinder::getPathNodes()
{
    return gridPath;
}

void gridPathFinder::resetPath()
{
    gridPath.clear();
}

// The state is a vector of {x, y, z, vx, vy, vz}
void gridPathFinder::stateTransit1(Eigen::VectorXd init_state, Eigen::VectorXd& final_state, Eigen::Vector3d um,
                                   double tau)
{
    // build the state transition matrix phi(tau)
    for (int i = 0; i < 3; ++i) phi(i, i + 3) = tau;

    // cout << "phi:\n" << phi << endl;

    // The build the integral terms
    Eigen::VectorXd integral(6);
    integral.head(3) = 0.5 * pow(tau, 2) * um;
    integral.tail(3) = tau * um;
    // cout << "integral:\n" << integral << endl;

    final_state = phi * init_state + integral;
    // cout << "The final state is:\n" << final_state.transpose() << endl << endl;
}

void gridPathFinder::jerkInputStateTransit(Eigen::VectorXd init_state, Eigen::VectorXd& final_state, Eigen::Vector3d um,
                                           double tau)
{
    Eigen::Vector3d x0, v0, a0;
    x0 = init_state.head(3);
    v0 = init_state.segment(3, 3);
    a0 = init_state.tail(3);

    // add tau in transit matrix
    for (int i = 0; i < 6; ++i) phi(i, i + 3) = tau;
    for (int i = 0; i < 3; ++i) phi(i, i + 6) = 0.5 * tau * tau;
    // std::cout << "phi:\n" << phi << std::endl;

    // integral
    Eigen::VectorXd integral(9);
    integral.head(3) = (1 / 6.0) * pow(tau, 3) * um;
    integral.segment(3, 3) = 0.5 * pow(tau, 2) * um;
    integral.tail(3) = tau * um;
    // std::cout << "integral:\n" << integral << std::endl;

    // cout << "init:" << init_state << endl;
    final_state = phi * init_state + integral;
    // cout << "final: " << final_state.transpose() << endl;
}

void gridPathFinder::getNeighbor(GridNodePtr current_node, GridNodePtr end_node, Neighbors& neighbors, double& vis_time)
{
    // applying different um and tau, we get a series of final state from current
    // state
    GridNodePtr nptr;
    Eigen::VectorXd state = current_node->state;

    Eigen::Vector3d um;
    double res = 1 / 2.0, time_res = 1 / 1.0;
    int pri_num = 0;
    for (double ax = -max_acc; ax <= max_acc + 1e-3; ax += max_acc * res)
        for (double ay = -max_acc; ay <= max_acc + 1e-3; ay += max_acc * res)
            for (double az = -max_acc; az <= max_acc + 1e-3; az += max_acc * res)
            {
                um << ax, ay, az;
                for (double tau = max_tau * time_res; tau <= max_tau + 1e-3; tau += max_tau * time_res)
                {
                    // state transit
                    Eigen::VectorXd x1;

                    stateTransit1(state, x1, um, tau);
                    // jerkInputStateTransit(state, x1, um, tau);

                    // display the transit
                    // ros::Time t1 = ros::Time::now();
                    // static int id_temp = 0;
                    // vector<Eigen::Vector3d> path_temp;
                    // for (double dt = 0; dt < tau + 1e-3; dt += 0.02)
                    // {
                    //     Eigen::VectorXd xt;
                    //     // jerkInputStateTransit(state, xt, um, dt);
                    //     stateTransit1(state, xt, um, dt);
                    //     path_temp.push_back(xt.head(3));
                    // }
                    // displayPathWithColor(path_temp, 0.01, 1, id_temp);
                    // ++id_temp;
                    // ros::Time t2 = ros::Time::now();
                    // vis_time += (t2 - t1).toSec();

                    // cout << "state:" << x1.transpose() << endl;

                    // stay in local range
                    Eigen::Vector3i idx1 = coord2gridIndex(x1.head(3));
                    if (idx1(0) < 0 || idx1(0) >= GLX_SIZE || idx1(1) < 0 || idx1(1) >= GLY_SIZE || idx1(2) < 0 ||
                        idx1(2) >= GLZ_SIZE)
                    {
                        // cout << "not in range" << endl;
                        break;
                    }

                    // collision free
                    nptr = GridNodeMap[idx1(0)][idx1(1)][idx1(2)];
                    if (nptr->occupancy > 0.5)
                    {
                        // cout << "obstacle" << endl;
                        break;
                    }

                    // vel feasible
                    Eigen::Vector3d v1 = x1.segment(3, 3);
                    if (fabs(v1(0)) > max_vel || fabs(v1(1)) > max_vel || fabs(v1(2)) > max_vel)
                    {
                        // cout << "vel end not feasible" << endl;
                        break;
                    }

                    // check if it is neighbor
                    Eigen::Vector3i diff = idx1 - current_node->index;
                    if (diff.norm() == 0) continue;
                    // cout << "neighbor:" << diff.transpose() << endl;
                    // caluculate f_score
                    double optimal_time;
                    KinoState candidate;
                    candidate.edge_cost = (um.squaredNorm() + w_time) * tau;
                    candidate.heu = getKinoDynamicHeu(x1, end_node->state, optimal_time);
                    candidate.optimal_time = optimal_time;
                    candidate.state = x1;
                    candidate.input = um;
                    candidate.duration = tau;

                    KinoState exist_neighbor;
                    // if not found, insert it directly
                    if (!neighbors.find(diff, exist_neighbor))
                    {
                        neighbors.add(diff, candidate);
                        // cout << "add new " << endl;
                    }
                    else  // compare the edge_cost + heu_cost
                    {
                        bool replace =
                            (candidate.edge_cost + candidate.heu) < (exist_neighbor.edge_cost + exist_neighbor.heu);
                        if (replace)
                        {
                            neighbors.erase(diff);
                            neighbors.add(diff, candidate);
                            // cout << "replace" << endl;
                        }
                    }
                }
            }
}

Eigen::Vector3i gridPathFinder::stringToIndex(std::string str)
{
    std::stringstream ss;
    ss << str;
    int idx;
    ss >> idx;
    int first = idx / 100;
    int second = (idx - first * 100) / 10;
    int third = (idx - first * 100 - second * 10);
    // cout << first << "," << second << "," << third << endl;
    Eigen::Vector3i diff;
    diff(0) = first - 1;
    diff(1) = second - 1;
    diff(2) = third - 1;
    return diff;
}

double gridPathFinder::getKinoDynamicHeu(GridNodePtr node1, GridNodePtr node2)
{
    const Vector3d dp = node2->state.head(3) - node1->state.head(3);
    const Vector3d v0 = node1->state.segment(3, 3);
    const Vector3d v1 = node2->state.segment(3, 3);

    double c1 = -36 * dp.dot(dp);
    double c2 = 24 * (v0 + v1).dot(dp);
    double c3 = -4 * (v0.dot(v0) + v0.dot(v1) + v1.dot(v1));
    double c4 = 0;
    double c5 = w_time;

    std::vector<double> ts = quartic(c5, c4, c3, c2, c1);
    double v_max = 2.0;
    double t_bar = (node1->state.head(3) - node2->state.head(3)).lpNorm<Infinity>() / v_max;
    ts.push_back(t_bar);

    double cost = std::numeric_limits<double>::max();
    double t_d = t_bar;

    for (auto t : ts)
    {
        if (t < t_bar) continue;
        double c = -c1 / 3 / t / t / t - c2 / 2 / t / t - c3 / t + w_time * t;
        if (c < cost)
        {
            cost = c;
            t_d = t;
        }
    }

    // node1->optimal_t = t_d;
    node1->optimal_time = t_d;

    return 1.0 * (1 + tie_breaker) * cost;
}

double gridPathFinder::getKinoDynamicHeu(Eigen::VectorXd x1, Eigen::VectorXd x2, double& optimal_time)
{
    const Vector3d dp = x2.head(3) - x1.head(3);
    const Vector3d v0 = x1.segment(3, 3);
    const Vector3d v1 = x2.segment(3, 3);

    double c1 = -36 * dp.dot(dp);
    double c2 = 24 * (v0 + v1).dot(dp);
    double c3 = -4 * (v0.dot(v0) + v0.dot(v1) + v1.dot(v1));
    double c4 = 0;
    double c5 = w_time;

    std::vector<double> ts = quartic(c5, c4, c3, c2, c1);
    double v_max = 2.0;
    double t_bar = (x1.head(3) - x2.head(3)).lpNorm<Infinity>() / v_max;
    ts.push_back(t_bar);

    double cost = std::numeric_limits<double>::max();
    double t_d = t_bar;

    for (auto t : ts)
    {
        if (t < t_bar) continue;
        double c = -c1 / 3 / t / t / t - c2 / 2 / t / t - c3 / t + w_time * t;
        if (c < cost)
        {
            cost = c;
            t_d = t;
        }
    }

    // node1->optimal_t = t_d;
    optimal_time = t_d;

    return 1.0 * (1 + tie_breaker) * cost;
}

vector<double> gridPathFinder::cubic(double a, double b, double c, double d)
{
    vector<double> dts;

    double a2 = b / a;
    double a1 = c / a;
    double a0 = d / a;

    double Q = (3 * a1 - a2 * a2) / 9;
    double R = (9 * a1 * a2 - 27 * a0 - 2 * a2 * a2 * a2) / 54;
    double D = Q * Q * Q + R * R;
    if (D > 0)
    {
        double S = std::cbrt(R + sqrt(D));
        double T = std::cbrt(R - sqrt(D));
        dts.push_back(-a2 / 3 + (S + T));
        return dts;
    }
    else if (D == 0)
    {
        double S = std::cbrt(R);
        dts.push_back(-a2 / 3 + S + S);
        dts.push_back(-a2 / 3 - S);
        return dts;
    }
    else
    {
        double theta = acos(R / sqrt(-Q * Q * Q));
        dts.push_back(2 * sqrt(-Q) * cos(theta / 3) - a2 / 3);
        dts.push_back(2 * sqrt(-Q) * cos((theta + 2 * M_PI) / 3) - a2 / 3);
        dts.push_back(2 * sqrt(-Q) * cos((theta + 4 * M_PI) / 3) - a2 / 3);
        return dts;
    }
}

vector<double> gridPathFinder::quartic(double a, double b, double c, double d, double e)
{
    vector<double> dts;

    double a3 = b / a;
    double a2 = c / a;
    double a1 = d / a;
    double a0 = e / a;

    vector<double> ys = cubic(1, -a2, a1 * a3 - 4 * a0, 4 * a2 * a0 - a1 * a1 - a3 * a3 * a0);
    double y1 = ys.front();
    double r = a3 * a3 / 4 - a2 + y1;
    if (r < 0) return dts;

    double R = sqrt(r);
    double D, E;
    if (R != 0)
    {
        D = sqrt(0.75 * a3 * a3 - R * R - 2 * a2 + 0.25 * (4 * a3 * a2 - 8 * a1 - a3 * a3 * a3) / R);
        E = sqrt(0.75 * a3 * a3 - R * R - 2 * a2 - 0.25 * (4 * a3 * a2 - 8 * a1 - a3 * a3 * a3) / R);
    }
    else
    {
        D = sqrt(0.75 * a3 * a3 - 2 * a2 + 2 * sqrt(y1 * y1 - 4 * a0));
        E = sqrt(0.75 * a3 * a3 - 2 * a2 - 2 * sqrt(y1 * y1 - 4 * a0));
    }

    if (!std::isnan(D))
    {
        dts.push_back(-a3 / 4 + R / 2 + D / 2);
        dts.push_back(-a3 / 4 + R / 2 - D / 2);
    }
    if (!std::isnan(E))
    {
        dts.push_back(-a3 / 4 - R / 2 + E / 2);
        dts.push_back(-a3 / 4 - R / 2 - E / 2);
    }

    return dts;
}

bool gridPathFinder::shotHeu(GridNodePtr node1, GridNodePtr node2)
{
    // cnt_shot++;
    // if (cnt_shot < N_max_shot)
    //     return false;
    // else
    // {
    //     cnt_shot = 0;
    //     N_max_shot = ceil((node1->state.head(3) - node2->state.head(3)).norm()
    //     / dis_shot * N_max);
    // }
    // TODO: need to check the velocity and acceleration feasibility

    const Vector3d p0 = node1->state.head(3);
    const Vector3d dp = node2->state.head(3) - p0;
    const Vector3d v0 = node1->state.segment(3, 3);
    const Vector3d v1 = node2->state.segment(3, 3);
    const Vector3d dv = v1 - v0;
    double t_d = node1->optimal_time;
    MatrixXd coef(3, 4);
    //  ****** now check the feasibility of the optimal polynomial, by using t_d

    while (true)
    {
        Vector3d a = 1.0 / 6.0 * (-12.0 / (t_d * t_d * t_d) * (dp - v0 * t_d) + 6 / (t_d * t_d) * dv);
        Vector3d b = 0.5 * (6.0 / (t_d * t_d) * (dp - v0 * t_d) - 2 / t_d * dv);

        coef.col(3) = a;
        coef.col(2) = b;
        coef.col(1) = v0;
        coef.col(0) = p0;

        // *** the OPTIMAL polynomial is : 1/6 * alpha * t^3 + 1/2 * beta * t^2 + v0
        // * t + p0; denote as : a*t^3 + b*t^2
        // + v0*t + p0
        Vector3d coord, vel, acc;
        VectorXd poly1d, t, polyv, polya;
        Vector3i index;

        double t_delta = 0.05;

        Eigen::MatrixXd Tm(4, 4);
        Tm << 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 3, 0, 0, 0, 0;

        bool feasible = true;

        // forward checking of trajectory
        for (double time = t_delta; time <= t_d; time += t_delta)
        {
            t = VectorXd::Zero(4);
            for (int j = 0; j < 4; j++) t(j) = pow(time, j);

            for (int dim = 0; dim < 3; dim++)
            {
                poly1d = coef.row(dim);
                coord(dim) = poly1d.dot(t);
                vel(dim) = (Tm * poly1d).dot(t);
                acc(dim) = (Tm * Tm * poly1d).dot(t);
            }

            index = coord2gridIndex(coord);

            if (index(0) < 0 || index(0) >= GLX_SIZE || index(1) < 0 || index(1) >= GLY_SIZE || index(2) < 0 ||
                index(2) >= GLZ_SIZE)
                return false;

            if (GridNodeMap[index(0)][index(1)][index(2)]->occupancy > 0.5)  // collision
                return false;

            // check kinodynamic feasibility
            for (int dim = 0; dim < 3; ++dim)
            {
                if (fabs(vel(dim) > 2.5))
                {
                    cout << "vel infeasible:" << vel << endl;
                    t_d *= fabs(vel(dim)) / 2.5;
                    cout << "new td:" << t_d << endl;
                    feasible = false;
                    break;
                }
            }

            if (!feasible) break;
        }

        if (feasible) break;
    }

    coef_shot = coef;
    t_shot = t_d;
    is_shot_succ = true;
    cout << "shot success!" << endl;
    return true;
}

vector<Eigen::Vector3d> gridPathFinder::getKinoTraj(double resolution)
{
    vector<Vector3d> state_list;

    GridNodePtr ptr = terminate_ptr;
    if (ptr == NULL)  // no path found
        return state_list;

    VectorXd xt, xu;

    // ROS_WARN("[hybridAstarSearch] check point's index");
    while (ptr->cameFrom != NULL)
    {
        Vector3d u = ptr->input;
        double duration = ptr->duration;
        // cout << "index:" << ptr->index.transpose() << ", duration:" << duration
        // << endl;

        xt = ptr->cameFrom->state;
        for (double t = duration; t >= -1e-3; t -= resolution)
        {
            stateTransit1(xt, xu, u, t);
            // jerkInputStateTransit(xt, xu, u, t);
            // cout << "t:" << t << ", state:" << xu.head(3).transpose() << endl;
            state_list.push_back(xu.head(3));
        }

        ptr = ptr->cameFrom;
    }

    reverse(state_list.begin(), state_list.end());

    if (is_shot_succ)  // add shoting heuristic trajectory to the kino traj list
    {
        Vector3d coord;
        VectorXd poly1d, t;

        for (double time = resolution; time <= t_shot; time += resolution)
        {
            for (int dim = 0; dim < 3; dim++)
            {
                poly1d = coef_shot.row(dim);
                t = VectorXd::Zero(4);

                for (int j = 0; j < 4; j++) t(j) = pow(time, j);

                coord(dim) = poly1d.dot(t);
            }

            state_list.push_back(coord);
        }
    }

    return state_list;
}

// input: ts', N(sample num)
// output: ts,  mat of 4 x (K+1)*(N+1) and each row for t, x, y, z
Eigen::MatrixXd gridPathFinder::getSamples(double& ts, int& K, int N)
{
    Eigen::MatrixXd samples;

    GridNodePtr ptr = terminate_ptr;
    if (ptr == NULL)
    {
        cout << "no path found, return null sample" << endl;
        return samples;
    }

    // cal the accumulated time of the path
    double T = t_shot;

    while (ptr->cameFrom != NULL)
    {
        T += ptr->duration;
        ptr = ptr->cameFrom;
    }

    cout << "accumulated time:" << T << endl;

    // cal ts, tm
    K = floor(T / ts);
    ts = T / (K + 1);
    double tm = ts / N;

    cout << "K:" << K << ", N:" << N << ", ts:" << ts << endl;

    // get samples
    bool in_shot = true;
    int Nt = 0;
    Eigen::VectorXd st((N + 1) * (K + 1)), sx((N + 1) * (K + 1)), sy((N + 1) * (K + 1)), sz((N + 1) * (K + 1));
    double T_accumulate = T, t = t_shot;

    while (true)
    {
        if (in_shot)  // cal one shot coordinate
        {
            Vector3d coord;
            VectorXd poly1d, tv;
            for (int dim = 0; dim < 3; dim++)
            {
                poly1d = coef_shot.row(dim);
                tv = VectorXd::Zero(4);

                for (int j = 0; j < 4; j++) tv(j) = pow(t, j);

                coord(dim) = poly1d.dot(tv);
            }

            st(Nt) = T_accumulate;
            sx(Nt) = coord(0);
            sy(Nt) = coord(1);
            sz(Nt) = coord(2);
            ++Nt;

            // segment connecting point must be added twice
            if (Nt % (N + 1) == 0 && Nt != (K + 1) * (N + 1))
            {
                st(Nt) = T_accumulate;
                sx(Nt) = coord(0);
                sy(Nt) = coord(1);
                sz(Nt) = coord(2);
                ++Nt;
            }

            // move to next sample
            t -= tm;
            T_accumulate -= tm;

            if (t < -1e-5)  // outside the range of path segment
            {
                in_shot = false;
                ptr = terminate_ptr;
                t += ptr->duration;
            }
        }
        else  // cal coordinate of normal path
        {
            Eigen::VectorXd xu;
            Vector3d u = ptr->input;
            Eigen::VectorXd xt = ptr->cameFrom->state;
            stateTransit1(xt, xu, u, t);

            st(Nt) = T_accumulate;
            sx(Nt) = xu(0);
            sy(Nt) = xu(1);
            sz(Nt) = xu(2);
            ++Nt;

            // segment connecting point must be added twice
            if (Nt % (N + 1) == 0 && Nt != (K + 1) * (N + 1))
            {
                st(Nt) = T_accumulate;
                sx(Nt) = xu(0);
                sy(Nt) = xu(1);
                sz(Nt) = xu(2);
                ++Nt;
            }

            // move to next sample
            t -= tm;
            T_accumulate -= tm;

            if (t < -1e-5)  // outside the range of path segment
            {
                if (ptr->cameFrom->cameFrom == NULL)  // reach the first node, finish all samples
                {
                    samples.resize(4, (K + 1) * (N + 1));
                    samples.row(0) = st.reverse();
                    samples.row(1) = sx.reverse();
                    samples.row(2) = sy.reverse();
                    samples.row(3) = sz.reverse();

                    return samples;
                }
                else  // not reach the first node
                {
                    ptr = ptr->cameFrom;
                    t += ptr->duration;
                }
            }
        }
    }
}