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

    w_time = 1.5;

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

    Vector3i start_idx = coord2gridIndex(start_pt);
    Vector3i end_idx = coord2gridIndex(end_pt);

    GridNodePtr startPtr = GridNodeMap[start_idx(0)][start_idx(1)][start_idx(2)];
    GridNodePtr endPtr = GridNodeMap[end_idx(0)][end_idx(1)][end_idx(2)];

    startPtr->state.head(3) = start_pt;
    startPtr->state.tail(3) = start_vel;
    endPtr->state.head(3) = end_pt;
    endPtr->state.tail(3) = end_vel;

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

        if (abs(current_node->index(0) - endPtr->index(0)) <= 2 &&
            abs(current_node->index(1) - endPtr->index(1)) <= 2 && abs(current_node->index(2) - endPtr->index(2)) <= 2)
        {
            // shotHeu(current_node, endPtr)
            ROS_WARN("[hybrid Astar]Reach goal..");
            // cout << "goal coord: " << endl << current_node->real_coord << endl;
            cout << "total number of iteration used in hybrid Astar: " << num_iter << endl;
            ros::Time time_2 = ros::Time::now();
            ROS_WARN("Time consume in hybrid A star path finding is %f", (time_2 - time_1).toSec());
            gridPath = retrievePath(current_node);
            terminate_ptr = current_node;
            has_path = true;
            return;
        }
        openSet.erase(openSet.begin());
        current_node->id = -1;  // move current_node node from open set to closed set.
        expandedNodes.push_back(current_node);

        // get the state of current node
        Eigen::VectorXd state = current_node->state;

        // cout << "\ncurrent state:" << state.transpose() << ", index:" << current_node->index.transpose() << endl;

        // get neighbor of this node
        std::multimap<std::string, KinoState> neighbors;
        getNeighbor(current_node, endPtr, neighbors);

        // iterate the neighbors
        std::multimap<std::string, KinoState>::iterator iter;

        for (iter = neighbors.begin(); iter != neighbors.end(); ++iter)
        {
            // get the neighbor node
            Eigen::Vector3i neighbor_idx, diff;
            diff = stringToIndex(iter->first);

            // neighbor diff is 012, so -1 is needed here.
            neighbor_idx(0) = current_node->index(0) + diff(0);
            neighbor_idx(1) = current_node->index(1) + diff(1);
            neighbor_idx(2) = current_node->index(2) + diff(2);

            neighbor_ptr = GridNodeMap[neighbor_idx(0)][neighbor_idx(1)][neighbor_idx(2)];

            double edge_cost = iter->second.edge_cost;
            double heu = iter->second.heu;
            double optimal_time = iter->second.optimal_time;
            double duration = iter->second.duration;
            Eigen::Vector3d input = iter->second.input;
            Eigen::VectorXd x1 = iter->second.state;

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
                    neighbor_ptr->nodeMapIt = openSet.insert(
                        make_pair(neighbor_ptr->fScore, neighbor_ptr));  // put neighbor in open set and record it.
                }
            }
            else if (neighbor_ptr->id != -1)  // not in the same grid and not in close set
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
                else if (tentative_gScore <= neighbor_ptr->gScore)  // already in open set, need compare and update
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
    Eigen::Vector3d x0, v0;
    x0 = init_state.head(3);
    v0 = init_state.tail(3);

    // cout << "The initial state is:\n" << init_state.transpose() << endl;
    // cout << "um:" << um.transpose() << ", tau:" << tau << endl;

    // build the state transition matrix phi(tau)
    Eigen::MatrixXd phi(6, 6);
    phi.block(3, 0, 3, 3) = Eigen::MatrixXd::Zero(3, 3);
    phi.block(0, 0, 3, 3) = phi.block(3, 3, 3, 3) = Eigen::MatrixXd::Identity(3, 3);
    phi.block(0, 3, 3, 3) = tau * Eigen::MatrixXd::Identity(3, 3);

    // cout << "phi:\n" << phi << endl;

    // The build the integral terms
    Eigen::VectorXd integral(6);
    integral.head(3) = 0.5 * pow(tau, 2) * um;
    integral.tail(3) = tau * um;
    // cout << "integral:\n" << integral << endl;

    final_state = phi * init_state + integral;
    // cout << "The final state is:\n" << final_state.transpose() << endl << endl;
}

void gridPathFinder::getNeighbor(GridNodePtr current_node, GridNodePtr end_node,
                                 std::multimap<std::string, KinoState>& neighbors)
{
    // applying different um and tau, we get a series of final state from current
    // state
    GridNodePtr nptr;
    Eigen::VectorXd state = current_node->state;

    Eigen::Vector3d um;
    double max_acc = 3.0;
    double max_tau = 0.3;
    double res = 1 / 2.0, time_res = 1 / 2.0;
    int pri_num = 0;
    for (double ax = -max_acc; ax <= max_acc + 1e-3; ax += max_acc * res)
        for (double ay = -max_acc; ay <= max_acc + 1e-3; ay += max_acc * res)
            for (double az = -max_acc; az <= max_acc + 1e-3; az += max_acc * res)
            {
                um << ax, ay, az;
                for (double tau = max_tau * time_res; tau <= max_tau + 1e-3; tau += max_tau * time_res)
                {
                    // first we should check the feasibility of the primitive,
                    // if it is not feasible, we drop it immediately
                    Eigen::VectorXd x1;
                    stateTransit1(state, x1, um, tau);

                    // stay in local range
                    Eigen::Vector3i idx1 = coord2gridIndex(x1.head(3));
                    if (idx1(0) < 0 || idx1(0) >= GLX_SIZE || idx1(1) < 0 || idx1(1) >= GLY_SIZE || idx1(2) < 0 ||
                        idx1(2) >= GLZ_SIZE)
                        break;

                    // collision free
                    nptr = GridNodeMap[idx1(0)][idx1(1)][idx1(2)];
                    if (nptr->occupancy > 0.5) break;

                    // dynamics feasible
                    Eigen::Vector3d v1 = x1.tail(3);
                    if (fabs(v1(0)) > 2.5 || fabs(v1(1)) > 2.5 || fabs(v1(2)) > 2.5) break;

                    // check if it is neighbor by idx difference, note that diff should be 012
                    Eigen::Vector3i diff = idx1 - current_node->index;

                    // the idx diff should only be +-1
                    // bool is_neighbor = diff.norm() != 0 && abs(diff(0)) <= 1 && abs(diff(1)) <= 1 && abs(diff(2)) <=
                    // 1;
                    bool is_neighbor = diff.norm() != 0;
                    if (!is_neighbor) continue;

                    // it is a neighbor! Save it. If there already have a neighbor with
                    // the same id, save the one have lower cost
                    string index =
                        to_string(int(diff(0)) + 1) + to_string(int(diff(1)) + 1) + to_string(int(diff(2)) + 1);

                    std::multimap<std::string, KinoState>::iterator iter = neighbors.find(index);

                    // caluculate f_score
                    // double cost = (um.squaredNorm() + rho) * tau;
                    double optimal_time;
                    KinoState candidate;
                    candidate.edge_cost = (um.squaredNorm() + w_time) * tau;
                    candidate.heu = getKinoDynamicHeu(x1, end_node->state, optimal_time);
                    candidate.optimal_time = optimal_time;
                    candidate.state = x1;
                    candidate.input = um;
                    candidate.duration = tau;

                    // if not found, insert it directly
                    if (iter == neighbors.end())
                    {
                        neighbors.insert(make_pair(index, candidate));
                    }
                    else  // already have one, compare the edge_cost + heu_cost
                    {
                        bool replace =
                            (candidate.edge_cost + candidate.heu) < (iter->second.edge_cost + iter->second.heu);
                        if (replace)
                        {
                            neighbors.erase(iter);
                            neighbors.insert(make_pair(index, candidate));
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
    const Vector3d v0 = node1->state.tail(3);
    const Vector3d v1 = node2->state.tail(3);

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
    const Vector3d v0 = x1.tail(3);
    const Vector3d v1 = x2.tail(3);

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
    cnt_shot++;
    if (cnt_shot < N_max_shot)
        return false;
    else
    {
        cnt_shot = 0;
        N_max_shot = ceil((node1->state.head(3) - node2->state.head(3)).norm() / dis_shot * N_max);
    }

    const Vector3d p0 = node1->state.head(3);
    const Vector3d dp = node2->state.head(3) - p0;
    const Vector3d v0 = node1->state.tail(3);
    const Vector3d v1 = node2->state.tail(3);
    const Vector3d dv = v1 - v0;
    const double t_d = node1->optimal_time;
    //  ****** now check the feasibility of the optimal polynomial, by using t_d

    Vector3d a = 1.0 / 6.0 * (-12.0 / (t_d * t_d * t_d) * (dp - v0 * t_d) + 6 / (t_d * t_d) * dv);
    Vector3d b = 0.5 * (6.0 / (t_d * t_d) * (dp - v0 * t_d) - 2 / t_d * dv);
    MatrixXd coef(3, 4);

    coef.col(3) = a;
    coef.col(2) = b;
    coef.col(1) = v0;
    coef.col(0) = p0;

    // *** the OPTIMAL polynomial is : 1/6 * alpha * t^3 + 1/2 * beta * t^2 + v0 * t + p0; denote as : a*t^3 + b*t^2 +
    // v0*t + p0
    Vector3d coord;
    VectorXd poly1d, t;
    Vector3i index;

    double t_delta = 0.05;

    for (double time = t_delta; time <= t_d; time += t_delta)
    {
        t = VectorXd::Zero(4);
        for (int j = 0; j < 4; j++) t(j) = pow(time, j);

        for (int dim = 0; dim < 3; dim++)
        {
            poly1d = coef.row(dim);
            coord(dim) = poly1d.dot(t);
        }

        index = coord2gridIndex(coord);

        if (index(0) < 0 || index(0) >= GLX_SIZE || index(1) < 0 || index(1) >= GLY_SIZE || index(2) < 0 ||
            index(2) >= GLZ_SIZE)
            return false;

        if (GridNodeMap[index(0)][index(1)][index(2)]->occupancy > 0.5)  // collision
            return false;
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
        // cout << "index:" << ptr->index.transpose() << ", duration:" << duration << endl;

        xt = ptr->cameFrom->state;
        for (double t = duration; t >= -1e-3; t -= resolution)
        {
            stateTransit1(xt, xu, u, t);
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
