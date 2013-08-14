
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "states.h"

#include "glog/logging.h"
#include "gflags/gflags.h"
#include "ceres/ceres.h"
#include "ceres/rotation.h"
DEFINE_bool(use_quaternions, false, "If true, uses quaternions to represent "
            "rotations. If false, angle axis is used.");

using namespace cerise;


bool DataLine::load(const std::string & line) {
    double f=0.0;
    int n=0;
    const char * lp = line.c_str();
    std::vector<double> dline;
    while (sscanf(lp, " %le%n", &f, &n)==1) {
        lp += n;
        dline.push_back(f);
    }
    if (dline.size() < 19) {
        LOG(ERROR) << "Not enough data filed on input line\n";
        return false;
    }
    timestamp = dline[0] + dline[1];
    depth = dline[8];
    vel = dline[9];
    vpred = dline[17];
    has_vel = round(dline[10]) > 0.5;
    dive_status = round(dline[18]);
    a[0] = dline[2];
    a[1] = dline[3];
    a[2] = dline[4];
    m[0] = dline[5];
    m[1] = dline[6];
    m[2] = dline[7];
    // conversion from ENU to NED (easier to manage with the compass and
    // accelerometer data in this convention)
    xyz[0] = dline[12];
    xyz[1] = dline[11];
    xyz[2] = -dline[13];
    // minus sign to convert from ENU to NED
    rpy[0] = -dline[14]*180./M_PI;
    rpy[1] = -dline[15]*180./M_PI;
    rpy[2] = -dline[16]*180./M_PI;
    return true;
}

bool GPSLine::load(const std::string & line) {
    double f=0.0;
    int n=0;
    const char * lp = line.c_str();
    std::vector<double> dline;
    while (sscanf(lp, " %le%n", &f, &n)==1) {
        lp += n;
        dline.push_back(f);
    }
    if (dline.size() < 12) {
        LOG(ERROR) << "Not enough data filed on gps input line\n";
        return false;
    }
    timestamp = dline[0] + dline[1];
    latitude = dline[2];
    longitude = dline[3];
    index = round(dline[4]) - 1;
    error = dline[5];
    easting = dline[6];
    northing = dline[7];
    zone = round(dline[8]);
    B[0] = dline[9];
    B[1] = dline[10];
    B[2] = dline[11];
    return true;
}

OptimisedOrientation::OptimisedOrientation() {
    for (size_t i=0;i<5;i++) {
        state[i] = 0.0;
    }
    rotation = state;
    if (FLAGS_use_quaternions) {
        propulsion = state + 4;
    } else {
        propulsion = state + 3;
    }
}

OptimisedOrientation::OptimisedOrientation(const OptimisedOrientation & oo) {
    dl_index = oo.dl_index;
    for (size_t i=0;i<5;i++) {
        state[i] = oo.state[i];
    }
    rotation = state;
    if (FLAGS_use_quaternions) {
        propulsion = state + 4;
    } else {
        propulsion = state + 3;
    }
}

bool OptimisedOrientation::load(const std::string & line, bool use_quaternion) {
    double f=0.0;
    int n=0;
    const char * lp = line.c_str();
    std::vector<double> dline;
    while (sscanf(lp, " %le%n", &f, &n)==1) {
        lp += n;
        dline.push_back(f);
    }
    if (dline.size() < (5 + (use_quaternion?1:0))) {
        LOG(ERROR) << "Not enough data filed on orientation line\n";
        return false;
    }
    dl_index = (size_t)(round(dline[0]));
    if (FLAGS_use_quaternions) {
        if (use_quaternion) {
            for (size_t i=0;i<4;i++) rotation[i] = dline[i+1];
            propulsion[0] = dline[5];
        } else {
            double qlog[3];
            for (size_t i=0;i<3;i++) qlog[i] = dline[i+1];
            ceres::AngleAxisToQuaternion(qlog,rotation);
            propulsion[0] = dline[4];
        }
    } else {
        if (use_quaternion) {
            double quat[4];
            for (size_t i=0;i<4;i++) quat[i] = dline[i+1];
            ceres::QuaternionToAngleAxis(quat,rotation);
            propulsion[0] = dline[5];
        } else {
            for (size_t i=0;i<3;i++) rotation[i] = dline[i+1];
            propulsion[0] = dline[4];
        }
    }
    return true;
}

std::string OptimisedOrientation::save() const {
    char buffer[1024];
    if (FLAGS_use_quaternions) {
        snprintf(buffer,1023,"%d %e %e %e %e %e",(int)dl_index,
                rotation[0],rotation[1],rotation[2],rotation[3],propulsion[0]);
    } else {
        snprintf(buffer,1023,"%d %e %e %e %e",(int)dl_index,
                rotation[0],rotation[1],rotation[2],propulsion[0]);
    }
    return std::string(buffer);
}

OptimisedOrientationSequence::OptimisedOrientationSequence() {
    use_quaternions = FLAGS_use_quaternions;
    Bscale = common_parameters;
    Kdepth = common_parameters + 3;
};

OptimisedOrientationSequence::OptimisedOrientationSequence(const OptimisedOrientationSequence & oos) {
    use_quaternions = FLAGS_use_quaternions;
    Bscale = common_parameters;
    Kdepth = common_parameters + 3;
    for (size_t i=0;i<4;i++) {
        common_parameters[i] = oos.common_parameters[i];
    }
    input_file = oos.input_file;
    states = oos.states;
}

void OptimisedOrientationSequence::initialise(const std::string & source_file, const std::vector<DataLine> & lines) 
{
    input_file = source_file;
    use_quaternions = FLAGS_use_quaternions;
    Bscale[0] = 100.0;
    Bscale[1] = 100.0;
    Bscale[2] = 100.0;
    Kdepth[0] = 2.0/600.0;
    states.clear();
    for (size_t i=0;i<lines.size();i++) {
        const DataLine & dl(lines[i]);
        OptimisedOrientation oo;
        double mat[9] = {1,0,0,0,1,0,0,0,1};
        oo.dl_index = i;
        ceres::EulerAnglesToRotationMatrix<double>(dl.rpy,3,mat);
        if (FLAGS_use_quaternions) {
            double log[3];
            ceres::RotationMatrixToAngleAxis<double>(mat,log);
            ceres::AngleAxisToQuaternion<double>(log,oo.rotation);
        } else {
            ceres::RotationMatrixToAngleAxis<double>(mat,oo.rotation);
        }
        states.push_back(oo);
    }
}

bool OptimisedOrientationSequence::load(const std::string & filename)
{
    states.clear();
    FILE * fp = fopen(filename.c_str(),"r");
    while (!feof(fp)) {
        char line[4096] = {0,};
        if (fgets(line,4095, fp)!= NULL) {
            std::string l(line);
            if (l.substr(0,5) == "#Inpu") {
                char iname[1024];
                sscanf(line,"#Input %s",iname);
                input_file = iname;
            } else if (l.substr(0,5) == "#UseQ") {
                int use_q = 0;
                sscanf(line,"#UseQuaternion %d",&use_q);
                use_quaternions = use_q;
            } else if (l.substr(0,5) == "#Bsca") {
                sscanf(line,"#Bscale %le %le %le",Bscale+0,Bscale+1,Bscale+2);
            } else if (l.substr(0,5) == "#Kdep") {
                sscanf(line,"#Kdepth %le",Kdepth);
            } else if (line[0] == '#') {
                continue;
            } else {
                OptimisedOrientation oo;
                if (oo.load(l,use_quaternions)) {
                    states.push_back(oo);
                }
            }
        }
    }
    use_quaternions = FLAGS_use_quaternions;
    return true;
}

bool OptimisedOrientationSequence::save(const std::string & filename) const 
{
    FILE * fp = fopen(filename.c_str(),"w");
    if (!fp) {
        LOG(ERROR) << "Can't write file: '"<<filename<<"'\n";
        return false;
    }
    fprintf(fp,"#Input %s\n",input_file.c_str());
    fprintf(fp,"#UseQuaternion %d\n",use_quaternions);
    fprintf(fp,"#Bscale %e %e %e\n",Bscale[0],Bscale[1],Bscale[2]);
    fprintf(fp,"#Kdepth %e\n",Kdepth[0]);
    fprintf(fp,"#Header index  R0 R1 R2 [R3] Prop\n");
    for (size_t i=0;i<states.size();i++) {
        fprintf(fp,"%s\n",states[i].save().c_str());
    }
    fclose(fp);
    return true;
}


OptimisedPosition::OptimisedPosition() {
    for (size_t i=0;i<5;i++) {
        state[i] = 0.0;
    }
    position = state;
    velocity = state+3;
    dzdt = state+4;
}

OptimisedPosition::OptimisedPosition(const OptimisedPosition & op) {
    dl_index = op.dl_index;
    for (size_t i=0;i<5;i++) {
        state[i] = op.state[i];
    }
    position = state;
    velocity = state+3;
    dzdt = state+4;
}

bool OptimisedPosition::load(const std::string & line) {
    double f=0.0;
    int n=0;
    const char * lp = line.c_str();
    std::vector<double> dline;
    while (sscanf(lp, " %le%n", &f, &n)==1) {
        lp += n;
        dline.push_back(f);
    }
    if (dline.size() < 5) {
        LOG(ERROR) << "Not enough data filed on Position line\n";
        return false;
    }
    dl_index = (size_t)(round(dline[0]));
    for (size_t i=0;i<5;i++) {
        state[i] = dline[i+1];
    }
    return true;
}

std::string OptimisedPosition::save() const {
    char buffer[1024];
    snprintf(buffer,1023,"%d %e %e %e %e %e",(int)dl_index,
            state[0],state[1],state[2],state[3],state[4]);
    return std::string(buffer);
}

OptimisedPositionSequence::OptimisedPositionSequence() {
    use_quaternions = FLAGS_use_quaternions;
    current = common_parameters;
    sensor_rotation = common_parameters + 2;
};

OptimisedPositionSequence::OptimisedPositionSequence(const OptimisedPositionSequence & ops) {
    use_quaternions = FLAGS_use_quaternions;
    current = common_parameters;
    sensor_rotation = common_parameters + 2;
    for (size_t i=0;i<6;i++) {
        common_parameters[i] = ops.common_parameters[i];
    }
    input_file = ops.input_file;
    states = ops.states;
}

void OptimisedPositionSequence::initialise(const std::string & source_file, const std::vector<DataLine> & lines) 
{
    input_file = source_file;
    use_quaternions = FLAGS_use_quaternions;
    current[0] = 0.0;
    current[1] = 0.0;
    if (use_quaternions) {
        double log[3] = {0,0,0};
        ceres::AngleAxisToQuaternion(log,sensor_rotation);
    } else {
        sensor_rotation[0] = sensor_rotation[1] = sensor_rotation[2] = 0.0;
    }
    states.clear();
    for (size_t i=0;i<lines.size();i++) {
        const DataLine & dl(lines[i]);
        OptimisedPosition op;
        op.dl_index = i;
        op.position[0] = dl.xyz[0];
        op.position[1] = dl.xyz[1];
        op.position[2] = dl.xyz[2];
        op.velocity[0] = dl.vpred; // predicted in the matlab script
        if (i > 1) {
            double dt = (lines[i].timestamp - lines[i-1].timestamp) * 24 * 3600;
            assert(dt > 0);
            op.dzdt[0] = -(lines[i].depth - lines[i-1].depth)/dt;
        } else {
            op.dzdt[0] = 0;
        }
        states.push_back(op);
    }
}

bool OptimisedPositionSequence::load(const std::string & filename)
{
    states.clear();
    FILE * fp = fopen(filename.c_str(),"r");
    while (!feof(fp)) {
        char line[4096] = {0,};
        if (fgets(line,4095, fp)!= NULL) {
            std::string l(line);
            if (l.substr(0,5) == "#Inpu") {
                char iname[1024];
                sscanf(line,"#Input %s",iname);
                input_file = iname;
            } else if (l.substr(0,5) == "#Orie") {
                char iname[1024];
                sscanf(line,"#Orientation %s",iname);
                orientation_file = iname;
            } else if (l.substr(0,5) == "#GPS ") {
                char iname[1024];
                sscanf(line,"#GPS %s",iname);
                gps_file = iname;
            } else if (l.substr(0,5) == "#UseQ") {
                int use_q = 0;
                sscanf(line,"#UseQuaternion %d",&use_q);
                use_quaternions = use_q;
            } else if (l.substr(0,5) == "#Curr") {
                sscanf(line,"#Current %le %le ",current+0,current+1);
            } else if (l.substr(0,5) == "#Rota") {
                double rot[4];
                if (use_quaternions) {
                    sscanf(line,"#Rotation %le %le %le %le",rot+0,rot+1,rot+2,rot+3);
                } else {
                    sscanf(line,"#Rotation %le %le %le",rot+0,rot+1,rot+2);
                }
                if (FLAGS_use_quaternions && !use_quaternions) {
                    ceres::AngleAxisToQuaternion(rot,sensor_rotation);
                } else if (!FLAGS_use_quaternions && use_quaternions) {
                    ceres::QuaternionToAngleAxis(rot,sensor_rotation);
                } else {
                    for (size_t i=0;i<4;i++) sensor_rotation[i] = rot[i];
                }
            } else if (line[0] == '#') {
                continue;
            } else {
                OptimisedPosition oo;
                if (oo.load(l)) {
                    states.push_back(oo);
                }
            }
        }
    }
    use_quaternions = FLAGS_use_quaternions;
    return true;
}

bool OptimisedPositionSequence::save(const std::string & filename) const 
{
    FILE * fp = fopen(filename.c_str(),"w");
    if (!fp) {
        LOG(ERROR) << "Can't write file: '"<<filename<<"'\n";
        return false;
    }
    fprintf(fp,"#Input %s\n",input_file.c_str());
    fprintf(fp,"#Orientation %s\n",orientation_file.c_str());
    fprintf(fp,"#GPS %s\n",gps_file.c_str());
    fprintf(fp,"#UseQuaternion %d\n",use_quaternions);
    fprintf(fp,"#Current %e %e \n",current[0],current[1]);
    if (use_quaternions) {
        fprintf(fp,"#Rotation %e %e %e %e\n",
                sensor_rotation[0],sensor_rotation[1],
                sensor_rotation[2],sensor_rotation[3]);
    } else {
        fprintf(fp,"#Rotation %e %e %e\n",
                sensor_rotation[0],sensor_rotation[1],sensor_rotation[2]);
    }
    fprintf(fp,"#Header index  X Y Z V dZ/dT\n");
    for (size_t i=0;i<states.size();i++) {
        fprintf(fp,"%s\n",states[i].save().c_str());
    }
    fclose(fp);
    return true;
}

