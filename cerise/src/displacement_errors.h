#ifndef ES_ROTATION_ERRORS_H
#define ES_ROTATION_ERRORS_H

#include "ceres/ceres.h"
#include "ceres/rotation.h"

namespace cerise {


    struct KinematicConstraint {
        KinematicConstraint(double dt, bool use_quat, double vpred, double *rot, double weight=1)
            : use_quaternions(use_quat), dt(dt), vpred(vpred), weight(weight) 
        {
            for (size_t i=0;i<(use_quat?4:3);i++) {
                rotation[i] = rot[i];
            }
        }

        // movement_parameters: 4D [ axis-angle rotation: 3, instantaneous propulsion acceleration: 1]
        // common_parameters: 4D [ magnetometer scale: 3, buoyancy gain: 1]
        template <typename T>
            bool operator()(const T* const current,
                    const T* const sensor_rotation, // Axis Angle
                    const T* const state1, 
                    const T* const state2, 
                    T* residuals) const {
                T rot[4];
                for (size_t i=0;i<(use_quaternions?4:3);i++) {
                    rot[i] = T(rotation[i]);
                }
                // T vdt_body[3] = {state1[3] * T(dt), T(0), T(0)};
                T vdt_body[3] = {T(vpred * dt), T(0), T(0)};
                T vdt_world[3];
                T vdt_world_corrected[3];
                if (use_quaternions) {
                    ceres::QuaternionRotatePoint(rot, vdt_body, vdt_world);
                    ceres::QuaternionRotatePoint(sensor_rotation, vdt_world, vdt_world_corrected);
                } else {
                    ceres::AngleAxisRotatePoint(rot, vdt_body, vdt_world);
                    ceres::AngleAxisRotatePoint(sensor_rotation, vdt_world, vdt_world_corrected);
                }

                // The error is the difference between the predicted and a position.
                residuals[0] = T(weight)*(state1[0] + vdt_world_corrected[0] + current[0] - state2[0]);
                residuals[1] = T(weight)*(state1[1] + vdt_world_corrected[1] + current[1] - state2[1]);
                residuals[2] = T(weight)*(state1[2] + vdt_world_corrected[2] - state2[2]);

                // In parallel, do the prediction of the vertical displacement
                // residuals[3] = T(weight)*(state1[2] + state1[4]*dt - state2[2]);

                return true;
            }

        bool use_quaternions;
        double dt;
        double vpred;
        double weight;
        double rotation[4];
    };

    struct VelocityError {
        VelocityError(double vel, double weight)
            : vel(vel), weight(weight) 
        {
        }

        template <typename T>
            bool operator()(const T* const state,
                    T* residuals) const {
                // The error is the difference between the predicted and a position.
                residuals[0] = T(weight)*(state[3] - T(vel));

                return true;
            }

        double vel;
        double weight;
    };

    struct DepthError {
        DepthError(double depth, double weight)
            : depth(depth), weight(weight) 
        {
        }

        template <typename T>
            bool operator()(const T* const position,
                    T* residuals) const {
                // The error is the difference between the predicted and a position.
                residuals[0] = T(weight)*(position[2] - (-T(depth)));

                return true;
            }

        double depth;
        double weight;
    };


    struct GPSError {
        GPSError(double gps_northing, double gps_easting, double weight)
            : gps_northing(gps_northing), gps_easting(gps_easting), weight(weight) 
        {
        }

        template <typename T>
            bool operator()(const T* const position,
                    T* residuals) const {
                // The error is the difference between the predicted and a position.
                residuals[0] = T(weight)*(position[0] - T(gps_northing));
                residuals[1] = T(weight)*(position[1] - T(gps_easting));

                return true;
            }

        double gps_northing;
        double gps_easting;
        double weight;
    };


};


#endif // ES_ROTATION_ERRORS_H
