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
    // cell_x_size_ * ((double)x_index + 0.5), cell_y_size_ * ((double)y_index + 0.5), cell_z_size_ * ((double)z_index +
    // 0.5)

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
    vector<GridNodePtr> visited_nodes;
    for (int i = 0; i < GLX_SIZE; i++)
        for (int j = 0; j < GLY_SIZE; j++)
            for (int k = 0; k < GLZ_SIZE; k++)
            {
                if (GridNodeMap[i][j][k]->id != 0)
                    // if(GridNodeMap[i][j][k]->id == -1)
                    visited_nodes.push_back(GridNodeMap[i][j][k]);
            }

    ROS_WARN("visited_nodes size : %d", visited_nodes.size());
    return visited_nodes;
}

/*bool gridPathFinder::minClearance()
{
    neighborPtr->occupancy > 0.5
}
*/
void gridPathFinder::AstarSearch(Eigen::Vector3d start_pt, Eigen::Vector3d end_pt)
{
    cout << "[Hybrid A Star]: begin" << endl;
    ros::Time time_1 = ros::Time::now();
    GridNodePtr startPtr = pos2gridNodePtr(start_pt);
    GridNodePtr endPtr = pos2gridNodePtr(end_pt);

    openSet.clear();

    GridNodePtr neighbor_ptr = NULL;
    GridNodePtr current_node = NULL;

    // initialize start and end point
    // Eigen::Vector3i start_idx = coord2gridIndex(start_pt);
    // Eigen::Vector3i end_idx = coord2gridIndex(end_pt);
    // GridNodePtr startPtr, endPtr;
    // startPtr = GridNodeMap[start_idx[0]][start_idx[1]][start_idx[2]];
    // endPtr = GridNodeMap[end_idx[0]][end_idx[1]][end_idx[2]];

    startPtr->gScore = 0;
    startPtr->fScore = getHeu(startPtr, endPtr);
    startPtr->id = 1;  // put start node in open set
    startPtr->coord = start_pt;
    openSet.insert(make_pair(startPtr->fScore, startPtr));  // put start in open set

    double tentative_gScore;
    int num_iter = 0;

    while (!openSet.empty())
    {
        ++num_iter;
        current_node = openSet.begin()->second;

        if (current_node->index(0) == endPtr->index(0) && current_node->index(1) == endPtr->index(1) &&
            current_node->index(2) == endPtr->index(2))
        {
            ROS_WARN("[Astar]Reach goal..");
            // cout << "goal coord: " << endl << current_node->real_coord << endl;
            cout << "total number of iteration used in Astar: " << num_iter << endl;
            ros::Time time_2 = ros::Time::now();
            ROS_WARN("Time consume in A star path finding is %f", (time_2 - time_1).toSec());
            gridPath = retrievePath(current_node);
            return;
        }
        openSet.erase(openSet.begin());
        current_node->id = -1;  // move current_node node from open set to closed set.
        expandedNodes.push_back(current_node);

        // get the state of current node
        Eigen::VectorXd state = current_node->state;
        cout << "\ncurrent state:" << state.transpose() << ", index:" << current_node->index.transpose() << endl;

        // get neighbor of this node
        multimap<std::string, pair<double, Eigen::VectorXd>> neighbors;
        getNeighbor(current_node, neighbors);

        multimap<std::string, pair<double, Eigen::VectorXd>>::iterator iter;
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

            // check the occupancy
            if (neighbor_ptr->occupancy > 0.5) continue;

            // skip node in the close set
            if (neighbor_ptr->id == -1) continue;

            cout << "nei state:" << iter->second.second.transpose() << ", nei idx:" << neighbor_ptr->index.transpose()
                 << endl;

            double cost = iter->second.first;
            Eigen::VectorXd x1 = iter->second.second;

            // to balance the heuristic and g, we need a lamda
            static double lamda = 0.2;
            tentative_gScore = current_node->gScore + lamda * cost;
            if (neighbor_ptr->id != 1)
            {
                // discover a new node
                neighbor_ptr->id = 1;
                neighbor_ptr->state = x1;
                neighbor_ptr->cameFrom = current_node;
                neighbor_ptr->gScore = tentative_gScore;
                // neighbor_ptr->fScore = neighbor_ptr->gScore + getHeu(neighbor_ptr, endPtr);
                neighbor_ptr->fScore = neighbor_ptr->gScore + 0.0;
                neighbor_ptr->nodeMapIt = openSet.insert(
                    make_pair(neighbor_ptr->fScore, neighbor_ptr));  // put neighbor in open set and record it.
                continue;
            }
            else if (tentative_gScore <= neighbor_ptr->gScore)
            {  // in open set and need update
                neighbor_ptr->state = x1;
                neighbor_ptr->cameFrom = current_node;
                neighbor_ptr->gScore = tentative_gScore;
                neighbor_ptr->fScore = tentative_gScore + getHeu(neighbor_ptr, endPtr);
                openSet.erase(neighbor_ptr->nodeMapIt);
                neighbor_ptr->nodeMapIt = openSet.insert(
                    make_pair(neighbor_ptr->fScore, neighbor_ptr));  // put neighbor in open set and record it.
            }
        }
    }

    ros::Time time_2 = ros::Time::now();
    ROS_WARN("Time consume in A star path finding is %f", (time_2 - time_1).toSec());
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

void gridPathFinder::getNeighbor(GridNodePtr current_node, multimap<string, pair<double, Eigen::VectorXd>>& neighbors)
{
    // applying different um and tau, we get a series of final state from current state
    Eigen::VectorXd state = current_node->state;

    Eigen::Vector3d um;
    double max_acc = 3.0;
    double max_tau = 0.5;
    double res = 1 / 2.0, time_res = 1 / 5.0;
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
                    Eigen::Vector3d v1 = x1.tail(3);
                    if (fabs(v1(0)) > 2.5 || fabs(v1(1)) > 2.5 || fabs(v1(2)) > 2.5) continue;

                    // Feasible! then we need to check whether it stay in local range
                    Eigen::Vector3i idx1 = coord2gridIndex(x1.head(3));
                    if (idx1(0) < 0 || idx1(0) >= GLX_SIZE || idx1(1) < 0 || idx1(1) >= GLY_SIZE || idx1(2) < 0 ||
                        idx1(2) >= GLZ_SIZE)
                        continue;

                    // check if it is neighbor, note that diff should be 012
                    Eigen::Vector3i diff = idx1 - current_node->index;
                    bool is_neighbor = !diff.norm() == 0 && abs(diff(0)) <= 1 && abs(diff(1)) <= 1 && abs(diff(2)) <= 1;
                    if (!is_neighbor) continue;

                    // it is a neighbor! Save it. If there already have a neighbor with the same id, save the one have
                    // lower cost
                    string index =
                        to_string(int(diff(0)) + 1) + to_string(int(diff(1)) + 1) + to_string(int(diff(2)) + 1);
                    multimap<string, pair<double, Eigen::VectorXd>>::iterator iter = neighbors.find(index);

                    // caluculate cost
                    double rho = 0.5;
                    double cost = (um.squaredNorm() + rho) * tau;
                    // if not found, insert it directly
                    if (iter == neighbors.end())
                    {
                        neighbors.insert(make_pair(index, make_pair(cost, x1)));
                    }
                    else  // already have one, compare the cost
                    {
                        bool replace = cost < iter->second.first;
                        if (replace)
                        {
                            neighbors.erase(iter);
                            neighbors.insert(make_pair(index, make_pair(cost, x1)));
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