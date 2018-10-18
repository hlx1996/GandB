#include <ros/ros.h>
#include <visualization_msgs/Marker.h>

#include <Eigen/Eigen>

#include <string>
#include <map>
#include <sstream>

using namespace std;

ros::Publisher primitive_pub;
ros::Publisher neighbor_pub;

int idx = 0;

// The state is a vector of {x, y, z, vx, vy, vz, ax, ay, az}
void stateTransit2(Eigen::VectorXd init_state, Eigen::VectorXd& final_state, Eigen::Vector3d um, double tau)
{
    Eigen::Vector3d x0, v0, a0;
    x0 = init_state.head(3);
    v0 = init_state.block(3, 0, 3, 1);
    a0 = init_state.tail(3);

    // cout << "The initial state is:\n" << init_state.transpose() << endl;
    // cout << "um:" << um.transpose() << ", tau:" << tau << endl;

    // build the state transition matrix phi(tau)
    Eigen::MatrixXd phi(9, 9);
    phi.block(3, 0, 6, 6) = Eigen::MatrixXd::Zero(6, 6);
    phi.block(0, 0, 3, 3) = phi.block(6, 6, 3, 3) = Eigen::MatrixXd::Identity(3, 3);
    phi.block(0, 3, 3, 3) = phi.block(3, 6, 3, 3) = tau * Eigen::MatrixXd::Identity(3, 3);
    phi.block(0, 6, 3, 3) = Eigen::MatrixXd::Identity(3, 3) * 0.5 * tau * tau;

    // cout << "phi:\n" << phi << endl;

    // The build the integral terms
    Eigen::VectorXd integral(9);
    integral.head(3) = (1 / 3.0) * pow(tau, 3) * um;
    integral.segment(3, 3) = 0.5 * pow(tau, 2) * um;
    integral.tail(3) = tau * um;
    // cout << "integral:\n" << integral << endl;

    final_state = phi * init_state + integral;
    // cout << "The final state is:\n" << final_state.transpose() << endl << endl;
}

// The state is a vector of {x, y, z, vx, vy, vz}
void stateTransit1(Eigen::VectorXd init_state, Eigen::VectorXd& final_state, Eigen::Vector3d um, double tau)
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

void generatePrimitive(const Eigen::VectorXd state, vector<Eigen::VectorXd>& finals, int layer)
{
    // Given a initial state, we apply a range of um and tau and get all the motion primitive
    // We store all the final state of the primitives in a vector
    // We visualize all these primitive in rviz
    // At the beginning I don't care about the feasibility. So I just uniformly sample some um and tau
    finals.clear();
    Eigen::Vector3d um;
    double dum = 3.0;
    double tau = 0.3;
    for (int i = -1; i <= 1; ++i)
        for (int j = -1; j <= 1; ++j)
            for (int k = 0; k <= 0; ++k)
            {
                um << double(i) * dum, double(j) * dum, double(k) * dum;
                // first we should check the feasibility of the primitive,
                // if it is not feasible, we drop it immediately
                Eigen::VectorXd x1;
                stateTransit1(state, x1, um, tau);
                Eigen::Vector3d v1 = x1.tail(3);
                if (fabs(v1(0)) > 2.5 || fabs(v1(1)) > 2.5 || fabs(v1(2)) > 2.5) continue;

                // Feasible! so we store the final state of this primitive
                cout << "final state is:" << x1.transpose() << endl;
                finals.push_back(x1);

                // The primitive is visualized
                visualization_msgs::Marker p;
                p.header.frame_id = "world";
                p.header.stamp = ros::Time::now();
                p.id = idx;

                p.type = visualization_msgs::Marker::SPHERE_LIST;
                p.action = visualization_msgs::Marker::ADD;

                p.pose.orientation.x = 0.0;
                p.pose.orientation.y = 0.0;
                p.pose.orientation.z = 0.0;
                p.pose.orientation.w = 1.0;
                p.scale.x = p.scale.y = p.scale.z = 0.01;

                if (layer == 1)
                {
                    p.color.a = p.color.r = 1.0;
                    p.color.g = p.color.b = 0.0;
                }
                else if (layer == 2)
                {
                    p.color.a = p.color.b = 1.0;
                    p.color.g = p.color.r = 0.0;
                }
                else if (layer == 3)
                {
                    p.color.a = p.color.g = 1.0;
                    p.color.b = p.color.r = 0.0;
                }
                for (double l = 0.0; l < tau; l += tau / 100)
                {
                    stateTransit1(state, x1, um, l);
                    geometry_msgs::Point pt;
                    pt.x = x1.head(3)(0);
                    pt.y = x1.head(3)(1);
                    pt.z = x1.head(3)(2);
                    p.points.push_back(pt);
                }

                primitive_pub.publish(p);
                ros::Duration(0.001).sleep();
                ++idx;
                cout << "idx: " << idx << endl;
            }
}

void generateNeighbor(const Eigen::VectorXd state)
{
    multimap<string, pair<double, Eigen::VectorXd>> neighbors;
    // we first get the center idx of x0, so that later we can use it to check neighbor
    Eigen::Vector3i center;
    double resolution = 0.2;
    double inv_res = 1 / resolution;
    Eigen::Vector3d origin;
    origin << 0.0, 0.0, 0.0;

    center(0) = int((state(0) - origin(0)) * inv_res);
    center(1) = int((state(1) - origin(1)) * inv_res);
    center(2) = int((state(2) - origin(2)) * inv_res);

    // applying different um and tau, we get a series of final state
    Eigen::Vector3d um;
    double max_acc = 3.0;
    double max_tau = 0.3;
    double res = 1 / 2.0;
    int pri_num = 0;
    for (double ax = -max_acc; ax <= max_acc + 1e-3; ax += max_acc * res)
        for (double ay = -max_acc; ay <= max_acc + 1e-3; ay += max_acc * res)
            for (double az = -max_acc; az <= max_acc + 1e-3; az += max_acc * res)
            {
                um << ax, ay, az;
                for (double tau = max_tau * res; tau <= max_tau + 1e-3; tau += max_tau * res)
                {
                    // first we should check the feasibility of the primitive,
                    // if it is not feasible, we drop it immediately
                    Eigen::VectorXd x1;
                    stateTransit1(state, x1, um, tau);
                    Eigen::Vector3d v1 = x1.tail(3);
                    if (fabs(v1(0)) > 2.5 || fabs(v1(1)) > 2.5 || fabs(v1(2)) > 2.5) continue;

                    // Feasible! but we still need to check whether it stay in neighbor grid
                    cout << "pri num:" << pri_num << "um: " << um.transpose() << " , tau" << tau << endl;
                    cout << "final state is:" << x1.transpose() << endl << endl;

                    visualization_msgs::Marker p;
                    p.header.frame_id = "world";
                    p.header.stamp = ros::Time::now();
                    p.id = pri_num;
                    p.type = visualization_msgs::Marker::SPHERE;
                    p.action = visualization_msgs::Marker::ADD;
                    p.pose.orientation.x = 0.0;
                    p.pose.orientation.y = 0.0;
                    p.pose.orientation.z = 0.0;
                    p.pose.orientation.w = 1.0;
                    p.scale.x = p.scale.y = p.scale.z = 0.01;

                    p.pose.position.x = x1(0);
                    p.pose.position.y = x1(1);
                    p.pose.position.z = x1(2);
                    p.color.r = p.color.a = 1.0;
                    p.color.b = p.color.g = 0.0;

                    primitive_pub.publish(p);
                    ++pri_num;
                    ros::Duration(0.001).sleep();

                    // now do the neighbor checking
                    Eigen::Vector3i idx1, diff;
                    idx1(0) = int((x1(0) - origin(0)) * inv_res);
                    idx1(1) = int((x1(1) - origin(1)) * inv_res);
                    idx1(2) = int((x1(2) - origin(2)) * inv_res);

                    diff = idx1 - center;

                    bool is_neighbor = !diff.norm() == 0 && abs(diff(0)) <= 1 && abs(diff(1)) <= 1 && abs(diff(2)) <= 1;

                    if (is_neighbor)
                    {
                        // calculate the cost of this neighbor primitive
                        double rho = 0.5;
                        double cost = (um.squaredNorm() + rho) * tau;

                        string index = to_string(int(diff(0))) + to_string(int(diff(1))) + to_string(int(diff(2)));
                        multimap<string, pair<double, Eigen::VectorXd>>::iterator iter = neighbors.find(index);
                        if (iter == neighbors.end())
                        {
                            neighbors.insert(make_pair(index, make_pair(cost, x1)));
                        }
                        else  // already have one, compare the cost
                        {
                            if (cost < iter->second.first)
                            {
                                neighbors.erase(iter);
                                neighbors.insert(make_pair(index, make_pair(cost, x1)));
                            }
                        }
                    }
                }
            }

    // visualize those neighbor points
    multimap<string, pair<double, Eigen::VectorXd>>::iterator iter;
    for (iter = neighbors.begin(); iter != neighbors.end(); ++iter)
    {
        Eigen::VectorXd x1 = iter->second.second.head(3);
        cout << "nei:" << x1.transpose() << endl;

        visualization_msgs::Marker p;
        p.header.frame_id = "world";
        p.header.stamp = ros::Time::now();
        p.id = pri_num;
        p.type = visualization_msgs::Marker::SPHERE;
        p.action = visualization_msgs::Marker::ADD;
        p.pose.orientation.x = 0.0;
        p.pose.orientation.y = 0.0;
        p.pose.orientation.z = 0.0;
        p.pose.orientation.w = 1.0;
        p.scale.x = p.scale.y = p.scale.z = 0.05;

        p.pose.position.x = x1(0);
        p.pose.position.y = x1(1);
        p.pose.position.z = x1(2);
        p.color.g = p.color.a = 1.0;
        p.color.b = p.color.r = 0.0;

        primitive_pub.publish(p);
        ++pri_num;
        ros::Duration(0.001).sleep();
    }
}

Eigen::Vector3i stringToIndex(string str)
{
    std::stringstream ss;
    ss << str;
    int idx;
    ss >> idx;
    int first = idx / 100;
    int second = (idx - first * 100) / 10;
    int third = (idx - first * 100 - second * 10);
    cout << first << "," << second << "," << third << endl;
}

int main(int argc, char** argv)
{
    ros::init(argc, argv, "astar");
    ros::NodeHandle node;

    primitive_pub = node.advertise<visualization_msgs::Marker>("primitive/path", 100);
    ros::Duration(0.5).sleep();

    cout << "test state transition" << endl;

    // given x0, v0, a0, we get the state by applying um for tau
    Eigen::VectorXd x0(6);
    x0 << 0.0, 0.0, 0.0, 1.0, 1.0, 0.0;

    // Eigen::Vector3d um;
    // um << 1.0, 1.0, 1.0;  // specify the input and duration
    // double tau = 1.0;

    // // calculate the final state
    // Eigen::VectorXd x1;
    // stateTransit1(x0, x1, um, tau);

    // // generate primitive
    // vector<Eigen::VectorXd> finals;
    // generatePrimitive(x0, finals, 1);

    // for (int i = 0; i < int(finals.size()); ++i)
    // {
    //     vector<Eigen::VectorXd> finals2;
    //     cout << i << " primitive:" << finals[i].transpose() << endl;
    //     generatePrimitive(finals[i], finals2, 2);
    // }

    generateNeighbor(x0);

    // ros::Duration(10.0).sleep();
    // test multimap here
    // multimap<Eigen::Vector3i, double> testmap;

    // Eigen::Vector3i vec1, vec2;
    // vec1 << 1, 1, 1;
    // vec2 << 1, 1, 1;

    // testmap.insert(make_pair(vec1, 2.0));
    // testmap.insert(make_pair(vec2, 1.0));

    // multimap<Eigen::Vector3i, double>::iterator iter;
    // for (iter = testmap.begin(); iter != testmap.end(); ++iter)
    // {
    //     cout << iter->first << ":\t" << iter->second << endl;
    // }

    // iter = testmap.find(vec1);
    // if (iter == testmap.end())
    //     cout << "not found" << endl;
    // else
    //     cout << iter->first << ":\t" << iter->second << endl;
    // string str = "-1-2-3";
    // cout << str[0] << ", " << str[1] << endl;

    // cout << "equal:" << (vec1 == vec2) << endl;

    stringToIndex("111");
    stringToIndex("012");
    stringToIndex("102");
    stringToIndex("100");

    return 0;
}