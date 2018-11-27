#ifndef _SDF_MAP_H
#define _SDF_MAP_H

#include <Eigen/Eigen>
#include <iostream>
#include <visualization_msgs/Marker.h>

using namespace std;

class SDFMap
{
  private:
    // data are saved in vector
    std::vector<int> occupancy_buffer;  // 0 is free, 1 is occupied
    std::vector<double> distance_buffer;
    std::vector<double> tmp_buffer1, tmp_buffer2;

    // map property
    Eigen::Vector3d origin, map_size;
    Eigen::Vector3d min_range, max_range;  // map range in pos
    Eigen::Vector3i grid_size;             // map range in index
    double resolution;
    Eigen::Vector3i min_vec, max_vec;  // the min and max updated range, unit is 1
    double truncated_distance = 20.0;

    bool isInMap(Eigen::Vector3d pos);
    void posToIndex(Eigen::Vector3d pos, Eigen::Vector3i& id);
    void indexToPos(Eigen::Vector3i id, Eigen::Vector3d& pos);

    template <typename F_get_val, typename F_set_val>
    void fillESDF(F_get_val f_get_val, F_set_val f_set_val, int start, int end, int dim);

  public:
    SDFMap(Eigen::Vector3d origin, double resolution, Eigen::Vector3d map_size);
    ~SDFMap()
    {
    }

    // occupancy management
    void setOccupancy(Eigen::Vector3d pos, int occ = 1);
    int getOccupancy(Eigen::Vector3d pos);
    void getOccupancyMarker(visualization_msgs::Marker& m, int id, Eigen::Vector4d color);

    // distance field management
    double getDistance(Eigen::Vector3d pos);
    void setUpdateRange(Eigen::Vector3d min_pos, Eigen::Vector3d max_pos);
    void updateESDF3d();
    void getESDFMarker(vector<visualization_msgs::Marker>& markers, int id, Eigen::Vector3d color);
    double getMaxDistance();
};

SDFMap::SDFMap(Eigen::Vector3d origin, double resolution, Eigen::Vector3d map_size)
{
    this->origin = origin;
    this->resolution = resolution;
    this->map_size = map_size;
    for (int i = 0; i < 3; ++i) grid_size(i) = ceil(map_size(i) / resolution);
    // cout << "grid num:" << grid_size.transpose() << endl;
    min_range = origin;
    max_range = origin + map_size;
    min_vec = Eigen::Vector3i::Zero();
    max_vec = grid_size - Eigen::Vector3i::Ones();

    // initialize size of buffer
    occupancy_buffer.resize(grid_size(0) * grid_size(1) * grid_size(2));
    distance_buffer.resize(grid_size(0) * grid_size(1) * grid_size(2));
    tmp_buffer1.resize(grid_size(0) * grid_size(1) * grid_size(2));
    tmp_buffer2.resize(grid_size(0) * grid_size(1) * grid_size(2));

    fill(distance_buffer.begin(), distance_buffer.end(), 10000);

    // initial map is free
    for (int i = 0; i < int(occupancy_buffer.size()); ++i) occupancy_buffer[i] = 0.0;
}

bool SDFMap::isInMap(Eigen::Vector3d pos)
{
    if (pos(0) < min_range(0) || pos(1) < min_range(1) || pos(2) < min_range(2))
    {
        cout << "less than min range!" << endl;
        return false;
    }

    if (pos(0) > max_range(0) || pos(1) > max_range(1) || pos(2) > max_range(2))
    {
        cout << "larger than max range!" << endl;
        return false;
    }

    return true;
}

void SDFMap::posToIndex(Eigen::Vector3d pos, Eigen::Vector3i& id)
{
    for (int i = 0; i < 3; ++i) id(i) = floor((pos(i) - origin(i)) / resolution);
}

void SDFMap::indexToPos(Eigen::Vector3i id, Eigen::Vector3d& pos)
{
    for (int i = 0; i < 3; ++i) pos(i) = (id(i) + 0.5) * resolution + origin(i);
}

void SDFMap::setOccupancy(Eigen::Vector3d pos, int occ)
{
    if (occ != 1 && occ != 0)
    {
        cout << "occ value error!" << endl;
        return;
    }

    if (!isInMap(pos)) return;

    Eigen::Vector3i id;
    posToIndex(pos, id);
    // cout << "id:" << id.transpose() << endl;

    // (x, y, z) -> x*ny*nz + y*nz + z
    // cout << "..." << id(0) * grid_size(1) * grid_size(2) + id(1) * grid_size(2) + id(2) << endl;
    // cout << "..." << occupancy_buffer.size() << endl;
    occupancy_buffer[id(0) * grid_size(1) * grid_size(2) + id(1) * grid_size(2) + id(2)] = occ;
}

int SDFMap::getOccupancy(Eigen::Vector3d pos)
{
    if (!isInMap(pos)) return -1;

    Eigen::Vector3i id;
    posToIndex(pos, id);

    // (x, y, z) -> x*ny*nz + y*nz + z
    return occupancy_buffer[id(0) * grid_size(1) * grid_size(2) + id(1) * grid_size(2) + id(2)];
}

void SDFMap::getOccupancyMarker(visualization_msgs::Marker& m, int id, Eigen::Vector4d color)
{
    m.header.frame_id = "world";
    m.id = id;
    m.type = visualization_msgs::Marker::CUBE_LIST;
    m.action = visualization_msgs::Marker::MODIFY;
    m.scale.x = resolution * 0.9;
    m.scale.y = resolution * 0.9;
    m.scale.z = resolution * 0.9;
    m.color.a = color(3);
    m.color.r = color(0);
    m.color.g = color(1);
    m.color.b = color(2);

    // iterate the map
    for (int x = 0; x < grid_size(0); ++x)
        for (int y = 0; y < grid_size(1); ++y)
            for (int z = 0; z < grid_size(2); ++z)
            {
                if (1 != occupancy_buffer[x * grid_size(1) * grid_size(2) + y * grid_size(2) + z]) continue;

                Eigen::Vector3d pos;
                indexToPos(Eigen::Vector3i(x, y, z), pos);

                geometry_msgs::Point p;
                p.x = pos(0);
                p.y = pos(1);
                p.z = pos(2);
                m.points.push_back(p);
            }
}

double SDFMap::getDistance(Eigen::Vector3d pos)
{
    if (!isInMap(pos)) return -1;

    Eigen::Vector3i id;
    posToIndex(pos, id);

    // (x, y, z) -> x*ny*nz + y*nz + z
    return distance_buffer[id(0) * grid_size(1) * grid_size(2) + id(1) * grid_size(2) + id(2)];
}

void SDFMap::setUpdateRange(Eigen::Vector3d min_pos, Eigen::Vector3d max_pos)
{
    if (!isInMap(min_pos)) min_pos = min_range;
    if (!isInMap(max_pos)) max_pos = max_range;

    posToIndex(min_pos, min_vec);
    posToIndex(max_pos - Eigen::Vector3d(resolution / 2, resolution / 2, resolution / 2), max_vec);
    cout << "min:" << min_vec.transpose() << ", max:" << max_vec << endl;
}

template <typename F_get_val, typename F_set_val>
void SDFMap::fillESDF(F_get_val f_get_val, F_set_val f_set_val, int start, int end, int dim)
{
    int v[grid_size(dim)];
    double z[grid_size(dim) + 1];

    int k = start;
    v[start] = start;
    z[start] = -std::numeric_limits<double>::max();
    z[start + 1] = std::numeric_limits<double>::max();

    for (int q = start + 1; q <= end; q++)
    {
        k++;
        double s;

        do
        {
            k--;
            s = ((f_get_val(q) + q * q) - (f_get_val(v[k]) + v[k] * v[k])) / (2 * q - 2 * v[k]);
            // ROS_INFO_STREAM("k: " << k << " s: " <<  s  << " z[k] " << z[k] << " v[k] " << v[k]);

        } while (s <= z[k]);

        k++;

        v[k] = q;
        z[k] = s;
        z[k + 1] = std::numeric_limits<double>::max();
    }

    k = start;

    for (int q = start; q <= end; q++)
    {
        while (z[k + 1] < q) k++;
        double val = (q - v[k]) * (q - v[k]) + f_get_val(v[k]);
        //      if(val < std::numeric_limits<_Scalar>::max())
        //  ROS_INFO_STREAM("val: " << val << " q: " << q << " v[k] " << v[k]);
        // if(val > truncation_distance_*truncation_distance_) val = std::numeric_limits<_Scalar>::max();
        f_set_val(q, val);
    }
}

void SDFMap::updateESDF3d()
{
    for (int x = min_vec[0]; x <= max_vec[0]; x++)
    {
        for (int y = min_vec[1]; y <= max_vec[1]; y++)
        {
            fillESDF(
                [&](int z) {
                    return occupancy_buffer[x * grid_size(1) * grid_size(2) + y * grid_size(2) + z] == 1 ?
                               0 :
                               std::numeric_limits<double>::max();
                },
                [&](int z, double val) { tmp_buffer1[x * grid_size(1) * grid_size(2) + y * grid_size(2) + z] = val; },
                min_vec[2], max_vec[2], 2);
        }
    }

    for (int x = min_vec[0]; x <= max_vec[0]; x++)
    {
        for (int z = min_vec[2]; z <= max_vec[2]; z++)
        {
            fillESDF(
                [&](int y) {
                    // cout << "get xyz:" << x << ", " << y << ", " << z << endl;
                    return tmp_buffer1[x * grid_size(1) * grid_size(2) + y * grid_size(2) + z];
                },
                [&](int y, double val) {
                    // cout << "set xyz:" << x << ", " << y << ", " << z << endl;
                    // cout << "index:" << x * grid_size(1) * grid_size(2) + y * grid_size(2) + z << endl;
                    // cout << "buffer length:" << tmp_buffer2.size() << endl;
                    tmp_buffer2[x * grid_size(1) * grid_size(2) + y * grid_size(2) + z] = val;
                },
                min_vec[1], max_vec[1], 1);
        }
    }

    for (int y = min_vec[1]; y <= max_vec[1]; y++)
    {
        for (int z = min_vec[2]; z <= max_vec[2]; z++)
        {
            fillESDF([&](int x) { return tmp_buffer2[x * grid_size(1) * grid_size(2) + y * grid_size(2) + z]; },
                     [&](int x, double val) {
                         distance_buffer[x * grid_size(1) * grid_size(2) + y * grid_size(2) + z] =
                             min(resolution * std::sqrt(val),
                                 distance_buffer[x * grid_size(1) * grid_size(2) + y * grid_size(2) + z]);
                     },
                     min_vec[0], max_vec[0], 0);
        }
    }

    min_vec = Eigen::Vector3i::Zero();
    max_vec = grid_size - Eigen::Vector3i::Ones();
}

void SDFMap::getESDFMarker(vector<visualization_msgs::Marker>& markers, int id, Eigen::Vector3d color)
{
    double max_dist = getMaxDistance();

    // get marker in several distance level
    const int level = 20;

    for (int i = 0; i < level; ++i)
    {
        visualization_msgs::Marker m;
        m.header.frame_id = "world";
        m.id = i + level * id;
        m.type = visualization_msgs::Marker::CUBE_LIST;
        m.action = visualization_msgs::Marker::ADD;
        m.scale.x = resolution * 0.9;
        m.scale.y = resolution * 0.9;
        m.scale.z = resolution * 0.9;
        m.color.r = color(0);
        m.color.g = color(1);
        m.color.b = color(2);

        // transparency and distance conversion
        double min_a = 0.05, max_a = 0.25;
        double da = (max_a - min_a) / (level - 1);
        m.color.a = max_a - da * i;
        // cout << "alpha:" << m.color.a << endl;

        // distance level
        double delta_d = max_dist / level;
        double min_d = i * delta_d;
        double max_d = (i + 1) * delta_d;
        // cout << "distance level:" << min_d << ", " << max_d << endl;

        // iterate the map
        for (int x = 0; x < grid_size(0); ++x)
            for (int y = 0; y < grid_size(1); ++y)
                for (int z = 0; z < grid_size(2) - 15; ++z)
                {
                    bool in_range = distance_buffer[x * grid_size(1) * grid_size(2) + y * grid_size(2) + z] < max_d &&
                                    distance_buffer[x * grid_size(1) * grid_size(2) + y * grid_size(2) + z] > min_d;
                    if (!in_range) continue;

                    Eigen::Vector3d pos;
                    indexToPos(Eigen::Vector3i(x, y, z), pos);

                    geometry_msgs::Point p;
                    p.x = pos(0);
                    p.y = pos(1);
                    p.z = pos(2);
                    m.points.push_back(p);
                }
        markers.push_back(m);
    }
}

double SDFMap::getMaxDistance()
{
    // get the max distance
    double max_dist = -1;
    for (int i = 0; i < int(distance_buffer.size()); ++i)
    {
        if (distance_buffer[i] > max_dist) max_dist = distance_buffer[i];
    }
    // cout << "Max distance is:" << max_dist << endl;
    return max_dist;
}

#endif