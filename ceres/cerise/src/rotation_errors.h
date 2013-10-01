#ifndef ES_ROTATION_ERRORS_H
#define ES_ROTATION_ERRORS_H

#include "ceres/ceres.h"
#include "ceres/rotation.h"

namespace cerise {


    struct AccelerometerError {
        AccelerometerError(double depth, double a_x, double a_y, double a_z, double weight=1)
            : depth(depth), a_x(a_x), a_y(a_y), a_z(a_z), weight(weight) 
        {
        }

        // movement_parameters: 4D [ axis-angle rotation: 3, instantaneous propulsion acceleration: 1]
        // common_parameters: 4D [ magnetometer scale: 3, buoyancy gain: 1]
        template <typename T>
            bool operator()(const T* const common_parameters,
                    const T* const rotation,
                    const T* const propulsion,
                    T* residuals) const {
                // camera[0,1,2] are the angle-axis rotation.
                T corrected_body_accel[3];
                corrected_body_accel[0] = T(a_x) - propulsion[0];
                corrected_body_accel[1] = T(a_y);
                corrected_body_accel[2] = T(a_z);
                T p[3];
                ceres::AngleAxisRotatePoint(rotation, corrected_body_accel, p);

                // The error is the difference between the predicted and a position.
                residuals[0] = T(weight)*(p[0] - T(G[0]));
                residuals[1] = T(weight)*(p[1] - T(G[1]));
                residuals[2] = T(weight)*(p[2] - (T(G[2]) + common_parameters[3]*T(depth)));

                return true;
            }

        static double G[3];

        double depth;
        double a_x;
        double a_y;
        double a_z;
        double weight;
    };

    struct MagnetometerError {
        MagnetometerError(double m_x, double m_y, double m_z, double weight)
            : m_x(m_x), m_y(m_y), m_z(m_z), weight(weight) 
        {
            b_x = B[0]; b_y = B[1]; b_z = B[2];
        }

        MagnetometerError(double m_x, double m_y, double m_z, 
                double b_x, double b_y, double b_z, 
                double weight)
            : b_x(b_x), b_y(b_y), b_z(b_z), 
            m_x(m_x), m_y(m_y), m_z(m_z), 
            weight(weight) 
        {
        }

        // movement_parameters: 3D [ axis-angle rotation: 3]
        // common_parameters: 4D [ magnetometer scale: 3, buoyancy gain: 1]
        template <typename T>
            bool operator()(const T* const common_parameters,
                    const T* const movement_parameters,
                    T* residuals) const {
                // camera[0,1,2] are the angle-axis rotation.
                T corrected_body_mag[3];
                corrected_body_mag[0] = T(m_x) * common_parameters[0];
                corrected_body_mag[1] = T(m_y) * common_parameters[1];
                corrected_body_mag[2] = T(m_z) * common_parameters[2];
                T p[3];
                ceres::AngleAxisRotatePoint(movement_parameters, corrected_body_mag, p);

                // The error is the difference between the predicted and a position.
                residuals[0] = T(weight)*(p[0] - T(b_x));
                residuals[1] = T(weight)*(p[1] - T(b_y));
                residuals[2] = T(weight)*(p[2] - T(b_z));

                return true;
            }

        static double B[3];

        double scale;
        double b_x;
        double b_y;
        double b_z;
        double m_x;
        double m_y;
        double m_z;
        double weight;
    };

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
        AccelerometerErrorQuat(double depth, double a_x, double a_y, double a_z, double weight=1)
            : depth(depth), a_x(a_x), a_y(a_y), a_z(a_z), weight(weight) 
        {
        }

        // movement_parameters: 5D [ Quaternion: 4, instantaneous propulsion acceleration: 1]
        // common_parameters: 4D [ magnetometer scale: 3, buoyancy gain: 1]
        template <typename T>
            bool operator()(const T* const common_parameters,
                    const T* const quaternion,
                    const T* const propulsion,
                    T* residuals) const {
                // camera[0,1,2] are the angle-axis rotation.
                T corrected_body_accel[3];
                corrected_body_accel[0] = T(a_x) - propulsion[0];
                corrected_body_accel[1] = T(a_y);
                corrected_body_accel[2] = T(a_z);
                T p[3];
                ceres::QuaternionRotatePoint(quaternion, corrected_body_accel, p);

                // The error is the difference between the predicted and a position.
                residuals[0] = T(weight)*(p[0] - T(G[0]));
                residuals[1] = T(weight)*(p[1] - T(G[1]));
                residuals[2] = T(weight)*(p[2] - (T(G[2]) + common_parameters[3]*T(depth)));

                return true;
            }

        static double G[3];

        double depth;
        double a_x;
        double a_y;
        double a_z;
        double weight;
    };

    struct MagnetometerErrorQuat {
        MagnetometerErrorQuat(double m_x, double m_y, double m_z, double weight)
            : m_x(m_x), m_y(m_y), m_z(m_z), weight(weight) 
        {
            b_x = B[0]; b_y = B[1]; b_z = B[2];
        }

        MagnetometerErrorQuat(double m_x, double m_y, double m_z, 
                double b_x, double b_y, double b_z, 
                double weight)
            : b_x(b_x), b_y(b_y), b_z(b_z), 
            m_x(m_x), m_y(m_y), m_z(m_z), 
            weight(weight) 
        {
        }


        // movement_parameters: 5D [ Quaternion: 4]
        // common_parameters: 4D [ magnetometer scale: 3, buoyancy gain: 1]
        template <typename T>
            bool operator()(const T* const common_parameters,
                    const T* const quaternion,
                    T* residuals) const {
                // camera[0,1,2] are the angle-axis rotation.
                T corrected_body_mag[3];
                corrected_body_mag[0] = T(m_x) * common_parameters[0];
                corrected_body_mag[1] = T(m_y) * common_parameters[1];
                corrected_body_mag[2] = T(m_z) * common_parameters[2];
                T p[3];
                ceres::QuaternionRotatePoint(quaternion, corrected_body_mag, p);

                // The error is the difference between the predicted and a position.
                residuals[0] = T(weight)*(p[0] - T(b_x));
                residuals[1] = T(weight)*(p[1] - T(b_y));
                residuals[2] = T(weight)*(p[2] - T(b_z));

                return true;
            }

        static double B[3];

        double scale;
        double b_x;
        double b_y;
        double b_z;
        double m_x;
        double m_y;
        double m_z;
        double weight;
    };

};


#endif // ES_ROTATION_ERRORS_H
