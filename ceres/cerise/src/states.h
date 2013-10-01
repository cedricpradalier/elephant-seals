#ifndef ES_STATES_H
#define ES_STATES_H

#include <string>
#include <vector>

namespace cerise {

    struct DataLine {
        double timestamp;
        double depth;
        bool has_vel;
        int dive_status;
        size_t dive_id;
        double vel;
        double vpred; // Initial velocity guess
        double a[3];
        double m[3];
        double xyz[3]; // Incremental position
        double rpy[3]; // Euler angles in degree

        bool load(const std::string & line) ;
        bool load_raw(const std::string & line) ;
    };

    struct GPSLine {
        double timestamp;
        double latitude;
        double longitude;
        int index;
        double error;
        double northing;
        double easting;
        int zone; // UTM Zone
        double B[3]; // Magnetic field

        bool load(const std::string & line) ;
    };


    struct OptimisedOrientation {
        size_t dl_index; // link to dataline index
        double state[5];
        double *rotation;
        double *propulsion;
        OptimisedOrientation(); 
        OptimisedOrientation(const OptimisedOrientation & oo); 

        bool load(const std::string & line, bool use_quaternion);

        std::string save() const;
        
    };

    struct OptimisedPosition {
        size_t dl_index; // link to dataline index
        double state[5];
        double *position; // 3D X,Y,Z NED
        double *velocity; // 1D linear velocity
        double *dzdt;     // 1D vertical velocity
        OptimisedPosition(); 
        OptimisedPosition(const OptimisedPosition & op); 

        bool load(const std::string & line);

        std::string save() const;
        
    };

    struct OptimisedOrientationSequence {
        std::string input_file;
        bool use_quaternions;
        double common_parameters[4];
        double *Bscale; // 3D vector of magnetometer gains
        double *Kdepth; // Buoyancy factor: measure accel = G + Kdepth*depth*[0;0;1]
        std::vector<OptimisedOrientation> states;

        OptimisedOrientationSequence();
        OptimisedOrientationSequence(const OptimisedOrientationSequence & oos); 

        void initialise(const std::string & source_file, const std::vector<DataLine> & lines);

        bool load(const std::string & filename);

        bool save(const std::string & filename) const;
    };

    struct OptimisedPositionSequence {
        std::string input_file;
        std::string orientation_file;
        std::string gps_file;
        bool use_quaternions;
        double common_parameters[6];
        double *current; // 2D (vx,vy) current velocity (averaged over the trip)
        double *sensor_rotation; // Sensor attachment orientation
        std::vector<OptimisedPosition> states;

        OptimisedPositionSequence();
        OptimisedPositionSequence(const OptimisedPositionSequence & oos); 

        void initialise(const std::string & source_file, const std::vector<DataLine> & lines);

        bool load(const std::string & filename);

        bool save(const std::string & filename) const;
    };

};



#endif // ES_STATES_H
