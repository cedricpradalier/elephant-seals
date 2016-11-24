

#include <stdlib.h>
#include <stdio.h>
#include "cgnuplot/CGnuplot.h"
#include "ceres/ceres.h"
#include "gflags/gflags.h"
#include "glog/logging.h"
#include "displacement_errors.h"
#include "states.h"

DEFINE_string(input, "", "Input File name");
DEFINE_string(orientations, "", "Orientation File name");
DEFINE_string(gps, "", "GPS File name");

DEFINE_string(trust_region_strategy, "levenberg_marquardt",
              "Options are: levenberg_marquardt, dogleg.");
DEFINE_string(dogleg, "traditional_dogleg", "Options are: traditional_dogleg,"
              "subspace_dogleg.");

DEFINE_bool(inner_iterations, false, "Use inner iterations to non-linearly "
            "refine each successful trust region step.");

DEFINE_string(blocks_for_inner_iterations, "automatic", "Options are: "
            "automatic, cameras, points, cameras,points, points,cameras");

DEFINE_string(linear_solver, "sparse_normal_cholesky", "Options are: "
              "sparse_schur, dense_schur, iterative_schur, sparse_normal_cholesky, "
              "dense_qr, dense_normal_cholesky and cgnr.");

DEFINE_string(preconditioner, "jacobi", "Options are: "
              "identity, jacobi, schur_jacobi, cluster_jacobi, "
              "cluster_tridiagonal.");

DEFINE_string(sparse_linear_algebra_library, "suite_sparse",
              "Options are: suite_sparse and cx_sparse.");

DEFINE_string(ordering, "automatic", "Options are: automatic, user.");

DEFINE_bool(robustify, false, "Use a robust loss function.");

DEFINE_double(eta, 1e-2, "Default value for eta. Eta determines the "
             "accuracy of each linear solve of the truncated newton step. "
             "Changing this parameter can affect solve performance.");

DEFINE_bool(use_local_parameterization, false, "For quaternions, use a local "
            "parameterization.");

DEFINE_int32(num_threads, 1, "Number of threads.");
DEFINE_int32(num_iterations, 50, "Number of iterations.");
DEFINE_double(max_solver_time, 1e32, "Maximum solve time in seconds.");
DEFINE_bool(nonmonotonic_steps, false, "Trust region algorithm can use"
            " nonmonotic steps.");

DEFINE_string(solver_log, "", "File to record the solver execution to.");



using namespace ceres;

namespace cerise{ 
    class OptimiseDisplacement {
        protected:
            cgnuplot::CGnuplot G;
            class DisplayCallback: public ceres::IterationCallback { 
                protected:
                    OptimiseDisplacement & problem;
                public: 
                    DisplayCallback(OptimiseDisplacement & p) : problem(p) {}
                    virtual ceres::CallbackReturnType operator()(const 
                            ceres::IterationSummary& summary) { 
                        problem.reportProgress();
                        // getchar();
                        return ceres::SOLVER_CONTINUE;
                    } 
            };

            DisplayCallback display;
            std::vector<DataLine> lines; 
            std::vector<GPSLine> gpslines; 
            OptimisedOrientationSequence oos;
            OptimisedPositionSequence ops;


            Problem problem;
            void SetLinearSolver(Solver::Options* options) {
                CHECK(StringToLinearSolverType(FLAGS_linear_solver,
                            &options->linear_solver_type));
                CHECK(StringToPreconditionerType(FLAGS_preconditioner,
                            &options->preconditioner_type));
                CHECK(StringToSparseLinearAlgebraLibraryType(
                            FLAGS_sparse_linear_algebra_library,
                            &options->sparse_linear_algebra_library));
                options->num_linear_solver_threads = FLAGS_num_threads;
            }

            void SetMinimizerOptions(Solver::Options* options) {
                options->max_num_iterations = FLAGS_num_iterations;
                options->minimizer_progress_to_stdout = true;
                options->num_threads = FLAGS_num_threads;
                options->eta = FLAGS_eta;
                options->function_tolerance = 3e-4;
                options->max_solver_time_in_seconds = FLAGS_max_solver_time;
                options->use_nonmonotonic_steps = FLAGS_nonmonotonic_steps;
                options->callbacks.push_back(&display);
                options->update_state_every_iteration = true; 

                CHECK(StringToTrustRegionStrategyType(FLAGS_trust_region_strategy,
                            &options->trust_region_strategy_type));
                CHECK(StringToDoglegType(FLAGS_dogleg, &options->dogleg_type));
                options->use_inner_iterations = FLAGS_inner_iterations;
            }

        public:
            OptimiseDisplacement() : display(*this) {}

            void optimise() {
                Solver::Options options;
                SetMinimizerOptions(&options);
                SetLinearSolver(&options);
                Solver::Summary summary;
                Solve(options, &problem, &summary);
                std::cout << summary.FullReport() << "\n";
            }

            bool load(const char * filename, const char * orientations, const char *gps, int lineLimit=-1) {
                lines.clear();
                FILE * fp = fopen(filename,"r");
                if (!fp) {
                    LOG(ERROR) << "Couldn't load data file '" << filename << "'\n";
                    return false;
                }
                while (!feof(fp)) {
                    if ((lineLimit>0) && ((signed)lines.size() >= lineLimit)) {
                        break;
                    }
                    char line[4096] = {0,};
                    if (fgets(line,4095, fp)!= NULL) {
                        DataLine dl;
                        if (dl.load(line)) {
                            lines.push_back(dl);
                        }
                    }
                }
                fclose(fp);
                printf("Loaded %d lines\n",(int)lines.size());

                fp = fopen(gps,"r");
                if (!fp) {
                    LOG(ERROR) << "Couldn't load gps file '" << gps << "'\n";
                    return false;
                }
                while (!feof(fp)) {
                    char line[4096] = {0,};
                    if (fgets(line,4095, fp)!= NULL) {
                        GPSLine dl;
                        if (dl.load(line)) {
                            gpslines.push_back(dl);
                        }
                    }
                }
                fclose(fp);
                printf("Loaded %d gps lines\n",(int)gpslines.size());
                assert(gpslines.size() >= 2);
                if (!oos.load(orientations)) {
                    LOG(ERROR) << "Couldn't load orientation file '" << orientations << "'\n";
                    return false;
                }
                lines[0].xyz[0] = gpslines[0].northing;
                lines[0].xyz[1] = gpslines[0].easting;
                lines[0].xyz[2] = -lines[0].depth;
                for (size_t i=1;i<lines.size();i++) {
                    double dt = (lines[i].timestamp-lines[i-1].timestamp)*24*3600;
                    double vdt_body[3] = {lines[i-1].vpred*dt,0,0};
                    double vdt_world[3];
                    if (oos.use_quaternions) {
                        ceres::QuaternionRotatePoint(oos.states[i-1].rotation,vdt_body,vdt_world);
                    } else {
                        ceres::AngleAxisRotatePoint(oos.states[i-1].rotation,vdt_body,vdt_world);
                    }
                    lines[i].xyz[0] = lines[i-1].xyz[0] + vdt_world[0];
                    lines[i].xyz[1] = lines[i-1].xyz[1] + vdt_world[1];
                    lines[i].xyz[2] = -lines[i].depth;
                }


                ops.initialise(filename,lines);
                ops.gps_file = gps;
                ops.orientation_file = orientations;
                

                for (size_t i=0;i<lines.size();i++) {
                    DataLine & dl(lines[i]);
                    OptimisedPosition & op(ops.states[i]);

                    CostFunction *cost_function;

                    cost_function = new AutoDiffCostFunction<cerise::DepthError,1,3>(
                            new cerise::DepthError(dl.depth, 1e3));
                    problem.AddResidualBlock(cost_function,NULL,op.state);

                    // LossFunction* loss_function;
                    // if (dl.has_vel) {
                    //     loss_function = FLAGS_robustify ? new HuberLoss(1.0) : NULL;
                    //     cost_function = new AutoDiffCostFunction<cerise::VelocityError,1,5>(
                    //             new cerise::VelocityError(dl.vel, 1.0));
                    //     problem.AddResidualBlock(cost_function,loss_function,op.state);
                    // }

                    if (i > 1) {
                        DataLine & dlp(lines[i-1]);
                        OptimisedPosition & opp(ops.states[i-1]);
                        OptimisedOrientation & oop(oos.states[i-1]);
                        double dt = (dl.timestamp - dlp.timestamp)*24*3600;

                        // Kinematics
                        if (ops.use_quaternions) {
                            cost_function = new AutoDiffCostFunction<cerise::KinematicConstraint,3,2,4,3,3>(
                                    new cerise::KinematicConstraint(dt,ops.use_quaternions,dlp.vpred,oop.rotation, 1.0e1));
                        } else {
                            cost_function = new AutoDiffCostFunction<cerise::KinematicConstraint,3,2,3,3,3>(
                                    new cerise::KinematicConstraint(dt,ops.use_quaternions,dlp.vpred,oop.rotation, 1.0e1));
                        }
                        problem.AddResidualBlock(cost_function,NULL,ops.current,ops.sensor_rotation, 
                                opp.state,op.state);

                    }

                    for (size_t i=0;i<gpslines.size();i++) {
                        GPSLine & gpl(gpslines[i]);
                        if (gpl.index < 0) {
                            continue;
                        }
                        assert(gpl.index < (signed)ops.states.size());
                        OptimisedPosition & op(ops.states[gpl.index]);
                        CostFunction *cost_function = new AutoDiffCostFunction<cerise::GPSError,2,3>(
                                new cerise::GPSError(gpl.northing,gpl.easting, 1e2));
                        problem.AddResidualBlock(cost_function,NULL,op.position);
                    }
                }
                if (ops.use_quaternions && FLAGS_use_local_parameterization) {
                    LocalParameterization* quaternion_parameterization = NULL;
                    quaternion_parameterization = new QuaternionParameterization;
                    problem.SetParameterization(ops.sensor_rotation, quaternion_parameterization);
                }
                reportProgress();
                getchar();
                return true;
            }

#if 1

            virtual void reportProgress() {
                if (ops.use_quaternions) {
                    printf("Current %.2f %.2f -- Sensor %.2f %.2f %.2f %.2f\n",
                            ops.current[0],ops.current[1],
                            ops.sensor_rotation[0],ops.sensor_rotation[1], 
                            ops.sensor_rotation[2], ops.sensor_rotation[3]);
                } else {
                    printf("Current %.2f %.2f -- Sensor %.2f %.2f %.2f\n",
                            ops.current[0],ops.current[1],
                            ops.sensor_rotation[0],ops.sensor_rotation[1], ops.sensor_rotation[2]);
                }
                ops.save("X");
                G.plot("set terminal wxt 1");
                G.plot("plot \"X\" u 0:(-$4) w l, \"%s\" u 0:9 w l",ops.input_file.c_str());
                G.plot("set terminal wxt 0");
                G.plot("plot \"X\" u 2:3 w l, \"%s\" u (-$13+%e):($12+%e) w l, \"%s\" u 8:7 w p ps 3",
                        ops.input_file.c_str(),gpslines[0].northing,gpslines[0].easting,ops.gps_file.c_str());
            }
#endif
    };
};


int main(int argc, char *argv[])
{
    google::ParseCommandLineFlags(&argc, &argv, true);
    google::InitGoogleLogging(argv[0]);
    if (FLAGS_input.empty() || FLAGS_orientations.empty() || FLAGS_gps.empty()) {
        LOG(ERROR) << "Usage: optimise_displacement --input preload.txt --gps gps.txt --orientations orientations.txt\n"
            << "\twhere preload and gps have been prepared by the prepare.m matlab script\n" 
            << "\tand orientations.txt has been generated by the optimise_orientation on the same dataset\n";
        return 1;
    }



    cerise::OptimiseDisplacement problem;

    if (problem.load(FLAGS_input.c_str(),FLAGS_orientations.c_str(),FLAGS_gps.c_str())) {

        problem.optimise();
        getchar();
    }
}


