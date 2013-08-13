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
            bool operator()(const T* const movement_parameters,
                    const T* const common_parameters,
                    T* residuals) const {
                // camera[0,1,2] are the angle-axis rotation.
                T corrected_body_accel[3];
                corrected_body_accel[0] = T(a_x) - movement_parameters[3];
                corrected_body_accel[1] = T(a_y);
                corrected_body_accel[2] = T(a_z);
                T p[3];
                ceres::AngleAxisRotatePoint(movement_parameters, corrected_body_accel, p);

                // The error is the difference between the predicted and a position.
                const T & K_depth = common_parameters[3];
                residuals[0] = T(weight)*(p[0] - G[0]);
                residuals[1] = T(weight)*(p[1] - G[1]);
                residuals[2] = T(weight)*(p[2] - (G[2] + K_depth*depth));

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
        }

        // movement_parameters: 4D [ axis-angle rotation: 3, instantaneous propulsion acceleration: 1]
        // common_parameters: 4D [ magnetometer scale: 3, buoyancy gain: 1]
        template <typename T>
            bool operator()(const T* const movement_parameters,
                    const T* const common_parameters,
                    T* residuals) const {
                // camera[0,1,2] are the angle-axis rotation.
                T corrected_body_mag[3];
                corrected_body_mag[0] = T(m_x) * movement_parameters[0];
                corrected_body_mag[1] = T(m_y) * movement_parameters[1];
                corrected_body_mag[2] = T(m_z) * movement_parameters[2];
                T p[3];
                ceres::AngleAxisRotatePoint(movement_parameters, corrected_body_mag, p);

                // The error is the difference between the predicted and a position.
                residuals[0] = T(weight)*(p[0] - B[0]);
                residuals[1] = T(weight)*(p[1] - B[1]);
                residuals[2] = T(weight)*(p[2] - B[2]);

                return true;
            }

        static double B[3];

        double scale;
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

        // movement_parameters: 4D [ axis-angle rotation: 3, instantaneous propulsion acceleration: 1]
        // common_parameters: 4D [ magnetometer scale: 3, buoyancy gain: 1]
        template <typename T>
            bool operator()(const T* const movement_parameters_1,
                    const T* const movement_parameters_2,
                    T* residuals) const {

                residuals[0] = T(weight) * (movement_parameters_2[3] - movement_parameters_1[3]);
                return true;
            }

        double weight;
    };
};


#endif // ES_ROTATION_ERRORS_H
