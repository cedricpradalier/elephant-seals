#ifndef ESO_MAT_OPS_H
#define ESO_MAT_OPS_H

#include <math.h>
#include <Eigen/Core.h>

inline Eigen::Matrix3f cross_mat(const Eigen::Vector3f & V) {
    Eigen::Matrix3f R;
    R << 0, -V(2), V(1),
      V(2), 0, -V(0),
      -V(1), V(0), 0;
    return R;
}

inline Eigen::Matrix3f rotx_mat(float roll) {
    Eigen::Matrix3f R;
    R << 1, 0, 0,
      0, cos(roll), -sin(roll),
      0, sin(roll), cos(roll);
    return R;
}

inline Eigen::Matrix3f roty_mat(float pitch) {
    Eigen::Matrix3f R;
    R << cos(roll), 0, sin(roll),
      0, 1, 0,
      -sin(roll), 0, cos(roll);
    return R;
}

inline Eigen::Matrix3f rotz_mat(float yaw) {
    Eigen::Matrix3f R;
    R << cos(roll), -sin(roll), 0,
      sin(roll), cos(roll), 0,
      0, 0, 1;
    return R;
}

inline Eigen::Matrix3f rpy_mat(float roll, float pitch, float yaw) {
    return rotz_mat(yaw)*roty_mat(pitch)*rotx_mat(roll);
}

inline Eigen::Matrix3f exp_mat(const Eigen::Vector3f & x) {
    Eigen::Matrix3f R;
    Eigen::Matrix3f cross = cross_mat(x);
    float theta = x.norm();
    R.setIdentity();
    if (fabs(theta)>1e-4) {
        R += cross*sin(theta)/theta + cross*cross*(1-cos(theta))/(theta*theta);
    }
    return R;
}

inline Eigen::Vector3f log_mat(const Eigen::Matrix3f & R) {
    Eigen::Vector3f V;
    float theta = acos((R.trace()-1)/2.0);
    V << R(2,1)-R(1,2), R(0,2)-R(2,0), R(1,0)-R(0,1);
    V /= 2;
    if (theta > 1e-4) {
        V *= theta/sin(theta);
    }
    return V;
}

#endif // ESO_MAT_OPS_H
