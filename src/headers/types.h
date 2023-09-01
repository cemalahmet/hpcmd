#ifndef __TYPES_H
#define __TYPES_H
#include <Eigen/Dense>

using Positions_t = Eigen::Array3Xd;
using Velocities_t = Eigen::Array3Xd;
using Forces_t = Eigen::Array3Xd;

using Vector3_t = Eigen::Vector3d;
using VectorX_t = Eigen::VectorXd;

using Names_t = std::vector<std::string>;

#endif // __TYPES_H
