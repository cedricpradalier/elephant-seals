#ifndef ES_ROTATION_ERRORS_H
#define ES_ROTATION_ERRORS_H

#include "ceres/ceres.h"
#include "ceres/rotation.h"
#include "boost/shared_ptr.hpp"

namespace cerise {

    struct SmoothnessConstraint {
        SmoothnessConstraint(double weight)
            : weight(weight)
        {
        }

        // movement_parameters: 1D [ instantaneous propulsion acceleration: 1]
        template <typename T>
            bool operator()(const T* const movement_parameters_1,
                    const T* const movement_parameters_2,
                    T* residuals) const {

                residuals[0] = T(weight) * (movement_parameters_2[0] - movement_parameters_1[0]);
                return true;
            }

        double weight;
    };

    struct AccelerometerErrorQuat {
        AccelerometerErrorQuat(double a_x, double a_y, double a_z, double weight=1)
            : a_x(a_x), a_y(a_y), a_z(a_z), weight(weight) 
        {
        }

        // movement_parameters: 4D [ Quaternion: 4]
        // common_parameters: 4D [ magnetometer scale: 3, buoyancy gain: 1]
        template <typename T>
            bool operator()(const T* const common_parameters,
                    const T* const quaternion,
                    T* residuals) const {
                T corrected_body_accel[3];
                corrected_body_accel[0] = T(a_x) * common_parameters[0];
                corrected_body_accel[1] = T(a_y) * common_parameters[0];
                corrected_body_accel[2] = T(a_z) * common_parameters[0];
                T p[3];
                ceres::QuaternionRotatePoint(quaternion, corrected_body_accel, p);

                // The error is the difference between the predicted and a position.
                residuals[0] = T(weight)*(p[0] - T(G[0]));
                residuals[1] = T(weight)*(p[1] - T(G[1]));
                residuals[2] = T(weight)*(p[2] - T(G[2]));

                return true;
            }

        static double G[3];

        double a_x;
        double a_y;
        double a_z;
        double weight;
    };

    struct MagnetometerErrorQuat {
        MagnetometerErrorQuat(bool northern,double m_x, double m_y, double m_z, 
                boost::shared_ptr<double> M, double weight)
            : northern(northern), m_x(m_x), m_y(m_y), m_z(m_z), M(M), weight(weight) 
        {
        }

        // movement_parameters: 4D [ Quaternion: 4]
        // common_parameters: 4D [ magnetometer scale: 3, buoyancy gain: 1]
        template <typename T>
            bool operator()(const T* const common_parameters,
                    const T* const quaternion,
                    T* residuals) const {
                T corrected_body_mag[3];
                corrected_body_mag[0] = T((m_x - M.get()[0])/M.get()[3]);
                corrected_body_mag[1] = T((m_y - M.get()[1])/M.get()[3]);
                corrected_body_mag[2] = T((m_z - M.get()[2])/M.get()[3]);
                T p[3];
                ceres::QuaternionRotatePoint(quaternion, corrected_body_mag, p);

                // The error is the difference between the predicted and a position.
                residuals[0] = T(weight)*(p[0] - cos(common_parameters[1]));
                residuals[1] = T(weight)*(p[1]); // - T(0.0));
                residuals[2] = T(weight)*(p[2] - sin(common_parameters[1])); 

                return true;
            }

        bool northern;
        double m_x;
        double m_y;
        double m_z;
        boost::shared_ptr<double> M;
        double weight;
    };

    struct ContinuityQuat {
        ContinuityQuat(double weight)
            : weight(weight) 
        {
        }

        template <typename T>
            bool operator()(const T* const q1,
                    const T* const q2,
                    T* residuals) const {
                T q2_inv[4], q_err[4];
                q2_inv[0] = q2[0];
                q2_inv[1] = -q2[1];
                q2_inv[2] = -q2[2];
                q2_inv[3] = -q2[3];
                ceres::QuaternionProduct(q2_inv, q1, q_err);
                ceres::QuaternionToAngleAxis(q_err,residuals);

                // The error is the difference between the predicted and a position.
                residuals[0] *= T(weight);
                residuals[1] *= T(weight);
                residuals[2] *= T(weight);

                return true;
            }

        double weight;
    };

    struct SphereConstraint {
        SphereConstraint(double vx, double vy, double vz):
            vx(vx), vy(vy), vz(vz)
        {
        }

        template <typename T>
            bool operator()(const T* const params,
                    T* residuals) const {

                // The error is the difference between the predicted and a position.
                residuals[0] = params[3] - ((T(vx)-params[0])*(T(vx)-params[0]) 
                        +(T(vy)-params[1])*(T(vy)-params[1]) 
                        +(T(vz)-params[2])*(T(vz)-params[2])); 

                return true;
            }

        double vx,vy,vz;
    };
};

#endif // ES_ROTATION_ERRORS_H
