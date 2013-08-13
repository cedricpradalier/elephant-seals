
#include "ros/ros.h" 
#include "tf/LinearMath/Matrix3x3.h" 
#include "es_optimizer/OptimisationProblem.h"
#include "cgnuplot/CGnuplot.h"

struct ORState : public StateDescription {

    ORState() {}
    virtual ~ORState() {}

    virtual size_t getDimension() {
        return 4;
    }


    virtual void step(const Eigen::VectorXf & delta) {
        R = exp_mat(delta.block(0,3)) * R;
        P += delta(3);
    }

    double timestamp;
    double depth;
    // Optimised variable
    double P;        // Propulsion acceleration
    Eigen::Matrix3f R; // Orientation
};

struct ORCState : public CommonState {

    ORCState() {}
    virtual ~ORCState() {}

    virtual size_t getDimension() {
        return 4;
    }

    virtual void step(const Eigen::VectorXf & delta) {
        Bscale += delta.block(0,3);
        k_depth += delta(3);
    }

    double Aref;
    double Bref;
    // Optimised variable
    double k_depth;         // Buoyancy correction gain to the accelerometer reading: A_buoy = [0;0;k_depth * depth]
    Eigen::Vector3f Bscale; // Magnetic field sensor gains
};

struct ORObservation : public ObservationDescription {
    ORObservation() {}
    virtual ~ORObservation() {}
    
    Eigen::Vector3f A;
    Eigen::Vector3f B;
};


struct MeasAccel : public ObservationError<ORCState,ORState,ORObservation> {
    Eigen::VectorXf error;
    Eigen::MatrixXf weights;
    Eigen::MatrixXf jacobian_state;
    Eigen::MatrixXf jacobian_common;

    MeasAccel(): error(3), weights(3,3), jacobian_state(3,4), jacobian_common(3,4) {
        weights.setIdentity();
        weights *= 10.0;
    }
    ~MeasAccel() {}

    virtual size_t getDimension() {return 3;}
    virtual const Eigen::VectorXf & getError() {return error;}
    virtual const Eigen::MatrixXf & getWeights() {return weights};
    virtual const Eigen::MatrixXf & getJacobian() {return jacobian_state};
    virtual const Eigen::MatrixXf & getJacobianWrtCommonState() {return jacobian_common};

    virtual void update(const Common & c, const State & s, const Observation & o) {
        Eigen::Vector3f propulsion; propulsion << s.P, 0, 0;
        Eigen::Vector3f buoyancy; buoyancy << 0, 0, s.depth * c.k_depth;
        error = s.R * (o.A - propulsion) - (c.Aref + buoyancy);
        jacobian_common(2,3) = -s.depth;
        jacobian_state.block(0,0,3,3) = - cross_mat(s.R * o.A);
        jacobian_state.block(0,3,3,1) = - s.R.block(0,0,3,1);
    }
}

struct MeasMag : public ObservationError<ORCState,ORState,ORObservation> {
    Eigen::VectorXf error;
    Eigen::MatrixXf weights;
    Eigen::MatrixXf jacobian_state;
    Eigen::MatrixXf jacobian_common;

    MeasMag(): error(3), weights(3,3), jacobian_state(3,4), jacobian_common(3,4) {
        weights.setIdentity();
        weights *= 0.1;
    }
    ~MeasMag() {}

    virtual size_t getDimension() {return 3;}
    virtual const Eigen::VectorXf & getError() {return error;}
    virtual const Eigen::MatrixXf & getWeights() {return weights};
    virtual const Eigen::MatrixXf & getJacobian() {return jacobian_state};
    virtual const Eigen::MatrixXf & getJacobianWrtCommonState() {return jacobian_common};

    virtual void update(const Common & c, const State & s, const Observation & o) {
        error = (c.Bscale.array() * o.B.array()).matrix() - s.R.transpose()*c.Bref;
        // jacobian_common(1:3,1:3) = diag(Bscale)
        jacobian_common(0,0) = c.Bscale(0);
        jacobian_common(1,1) = c.Bscale(1);
        jacobian_common(2,2) = c.Bscale(2);
        jacobian_state.block(0,0,3,3) = - s.R.transpose() * cross_mat(c.Bref);
    }
}

struct ContProp : public TransitionError<ORCState,ORState,ORObservation> {
    Eigen::VectorXf error;
    Eigen::MatrixXf weights;
    Eigen::MatrixXf jacobian_state;
    Eigen::MatrixXf jacobian_common;

    ContProp(): error(1), weights(1,1), jacobian_state(1,4), jacobian_common(1,4) {
        weights(0,0) = 1e3;
        jacobian_state(0,7) = 1;
        jacobian_state(0,3) = -1;
    }
    ~ContProp() {}

    virtual size_t getDimension() {return 3;}
    virtual const Eigen::VectorXf & getError() {return error;}
    virtual const Eigen::MatrixXf & getWeights() {return weights};
    virtual const Eigen::MatrixXf & getJacobian() {return jacobian_state};
    virtual const Eigen::MatrixXf & getJacobianWrtCommonState() {return jacobian_common};

    virtual void update(const Common & c, const State & s, const Observation & o) {
        error = s2.P - s1.P;
    }
}
typedef boost::shared_ptr<MeasAccel> MeasAccelPtr;
typedef boost::shared_ptr<MeasMag> MeasMagPtr;
typedef boost::shared_ptr<ContProp> ContPropPtr;

class OROptimisation : public OptimisationProblem<ORCState,ORState,ORObservation> {
    protected:
        CGnuplot G;

    public:
        OROptimisation() {
            observationErrors.push_back(MeasAccelPtr(new MeasAccel()));
            observationErrors.push_back(MeasMagPtr(new MeasMag()));
            transitionErrors.push_back(ContPropPtr(new ContProp()));

            // Extracted with igrf, from the matlab side
            common.Bref << 9.6920e+03, -1.4563e+04, -4.5355e+04;
            common.Aref << 0,0,-10;
        }
        ~OROptimisation() {}

        bool load(int lineLimit=-1) {
            common.k_depth = 2.0/600.0;
            common.Bscale << 100., 100., 100.;
            FILE * fp = fopen("data/preload_mat.txt","r");
            while (!feof(fp)) {
                char line[4096] = {0,};
                std::vector<double> dline;
                if (fgets(line,4095, fp)!= NULL) {
                    double f=0.0;
                    int n=0;
                    char * lp = line;
                    while (sscanf(lp, " %e%n", &f, &n)==1) {
                        lp += n;
                        dline.push_back(f);
                    }
                }
                ORState s; 
                s.timestamp = dline[0] + dline[1];
                s.depth = dline[9];
                s.P = 0;
                s.R.setIdentity();
                // alternatively, but does not seem necessary: s.R = rpy_mat(dline[15],dline[16],dline[17]);
                states.push_back(s);

                ORObservation o;
                o.A << dline[3], dline[4], dline[5];
                o.B << dline[6], dline[7], dline[8];
                if ((lineLimit>0) && ((signed)states.size() >= lineLimit)) {
                    break;
                }
            }
            fclose(fp);
            updateMatrices();
            printf("Loaded %d lines\n",states.size());
            return true;
        }

        virtual void reportProgress() {
            FILE * fa = fopen("Ap","w"), fb = fopen("Bp","w");
            for (size_t i=0;i<states.size();i++) {
                Eigen::Vector3f Ap, Bp, buoyancy;
                buoyancy << 0, 0, common.k_depth * states[i].depth;
                Ap = states[i].R * obs[i] - buoyancy;
                Bp = states[i].R * (common.Bscale.array() * obs[i].B.array()).matrix();
                double t = (states[i].timestamp - states[0].timestamp)*24*3600;
                fprintf(fa,"%e %e %e %e\n",t,Ap(0),Ap(1),Ap(2));
                fprintf(fb,"%e %e %e %e\n",t,Bp(0),Bp(1),Bp(2));
            }
            fclose(fa); fclose(fb);
            G.plot("set terminal x11 0");
            G.plot("plot \"Ap\" u 1:2 w l, \"Ap\" u 1:3 w l, \"Ap\" u 1:4 w l");
            G.plot("set terminal x11 1");
            G.plot("plot \"Bp\" u 1:2 w l, \"Bp\" u 1:3 w l, \"Bp\" u 1:4 w l");
        }
};


int main(int argc, char *argv[])
{

    OROptimisation problem;

    problem.load();
    
    problem.optimise();
}


