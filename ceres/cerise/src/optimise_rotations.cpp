

#include <stdlib.h>
#include <stdio.h>
#include "cgnuplot/CGnuplot.h"
#include "ceres/ceres.h"
#include "gflags/gflags.h"
#include "glog/logging.h"
#include "rotation_errors.h"

DEFINE_string(input, "", "Input File name");
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

DEFINE_bool(use_quaternions, false, "If true, uses quaternions to represent "
            "rotations. If false, angle axis is used.");

DEFINE_bool(use_local_parameterization, false, "For quaternions, use a local "
            "parameterization.");

DEFINE_bool(robustify, false, "Use a robust loss function.");

DEFINE_double(eta, 1e-2, "Default value for eta. Eta determines the "
             "accuracy of each linear solve of the truncated newton step. "
             "Changing this parameter can affect solve performance.");

DEFINE_int32(num_threads, 1, "Number of threads.");
DEFINE_int32(num_iterations, 50, "Number of iterations.");
DEFINE_double(max_solver_time, 1e32, "Maximum solve time in seconds.");
DEFINE_bool(nonmonotonic_steps, false, "Trust region algorithm can use"
            " nonmonotic steps.");

DEFINE_string(solver_log, "", "File to record the solver execution to.");



using namespace ceres;

class OptimiseRotation {
    protected:
        cgnuplot::CGnuplot G;
        struct DataLine {
            double timestamp;
            double depth;
            double a[3];
            double m[3];
            double rpy[3];

            double state[5];
        };

        class DisplayCallback: public ceres::IterationCallback { 
            protected:
                OptimiseRotation & problem;
            public: 
                DisplayCallback(OptimiseRotation & p) : problem(p) {}
                virtual ceres::CallbackReturnType operator()(const 
                        ceres::IterationSummary& summary) { 
                    problem.reportProgress();
                    // getchar();
                    return ceres::SOLVER_CONTINUE;
                } 
        };

        DisplayCallback display;
        double common_parameters[4];
        std::vector<DataLine> states; /* will be N*4 at the end */
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
            options->function_tolerance = 5e-4;
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
        OptimiseRotation() : display(*this) {}

        void optimise() {
            Solver::Options options;
            SetMinimizerOptions(&options);
            SetLinearSolver(&options);
            Solver::Summary summary;
            Solve(options, &problem, &summary);
            std::cout << summary.FullReport() << "\n";
        }

        bool load(const char * filename, int lineLimit=-1) {
            states.clear();
            /* Bscale */
            common_parameters[0] = 100.0;
            common_parameters[1] = 100.0;
            common_parameters[2] = 100.0;
            /*K_depth =*/ 
            common_parameters[3] = 2.0/600.0;
            FILE * fp = fopen(filename,"r");
            while (!feof(fp)) {
                char line[4096] = {0,};
                std::vector<double> dline;
                if (fgets(line,4095, fp)!= NULL) {
                    double f=0.0;
                    int n=0;
                    char * lp = line;
                    while (sscanf(lp, " %le%n", &f, &n)==1) {
                        lp += n;
                        dline.push_back(f);
                    }
                    assert(dline.size() >= 18);
                    DataLine dl;
                    dl.timestamp = dline[0] + dline[1];
                    dl.depth = dline[8];
                    dl.a[0] = dline[2];
                    dl.a[1] = dline[3];
                    dl.a[2] = dline[4];
                    dl.m[0] = dline[5];
                    dl.m[1] = dline[6];
                    dl.m[2] = dline[7];
                    dl.rpy[0] = -dline[14]*180./M_PI;
                    dl.rpy[1] = -dline[15]*180./M_PI;
                    dl.rpy[2] = -dline[16]*180./M_PI;
                    states.push_back(dl);
                    if ((lineLimit>0) && ((signed)states.size() >= lineLimit)) {
                        break;
                    }
                }
            }
            fclose(fp);
            printf("Loaded %d lines\n",(int)states.size());
                
            for (size_t i=0;i<states.size();i++) {
                DataLine & dl(states[i]);
                double mat[9] = {1,0,0,0,1,0,0,0,1};
                EulerAnglesToRotationMatrix<double>(dl.rpy,3,mat);
                if (FLAGS_use_quaternions) {
                    double log[3];
                    RotationMatrixToAngleAxis<double>(mat,log);
                    AngleAxisToQuaternion<double>(log,dl.state);
                } else {
                    RotationMatrixToAngleAxis<double>(mat,dl.state);
                }
                // printf("RPY %e %e %e LOG %e %e %e\n",
                //         dl.rpy[0],dl.rpy[1],dl.rpy[2],
                //         dl.state[0],dl.state[1],dl.state[2]);
                // If enabled use Huber's loss function.
                LossFunction* loss_function;
                CostFunction *cost_function;
                // Acceleration
                loss_function = FLAGS_robustify ? new HuberLoss(1.0) : NULL;
                if (FLAGS_use_quaternions) {
                    cost_function = new AutoDiffCostFunction<cerise::AccelerometerErrorQuat,3,4,4>(
                            new cerise::AccelerometerErrorQuat(dl.depth,dl.a[0],dl.a[1],dl.a[2], 10.0));
                } else {
                    cost_function = new AutoDiffCostFunction<cerise::AccelerometerError,3,4,4>(
                            new cerise::AccelerometerError(dl.depth,dl.a[0],dl.a[1],dl.a[2], 10.0));
                }
                problem.AddResidualBlock(cost_function,loss_function,common_parameters, states[i].state);

                // Magnetic field
                loss_function = FLAGS_robustify ? new HuberLoss(100.0) : NULL;
                if (FLAGS_use_quaternions) {
                    cost_function = new AutoDiffCostFunction<cerise::MagnetometerErrorQuat,3,4,4>(
                            new cerise::MagnetometerErrorQuat(dl.m[0],dl.m[1],dl.m[2], 0.01));
                } else {
                    cost_function = new AutoDiffCostFunction<cerise::MagnetometerError,3,4,4>(
                            new cerise::MagnetometerError(dl.m[0],dl.m[1],dl.m[2], 0.01));
                }
                problem.AddResidualBlock(cost_function,loss_function,common_parameters, states[i].state);

                if (i>0) {
                    // Smoothness constraint
                    if (FLAGS_use_quaternions) {
                        cost_function = new AutoDiffCostFunction<cerise::SmoothnessConstraintQuat,1,4,4>(
                                new cerise::SmoothnessConstraintQuat(1e3));
                    } else {
                        cost_function = new AutoDiffCostFunction<cerise::SmoothnessConstraint,1,4,4>(
                                new cerise::SmoothnessConstraint(1e3));
                    }
                    problem.AddResidualBlock(cost_function,NULL,states[i-1].state, states[i].state);
                }
                if (FLAGS_use_quaternions) {
                    LocalParameterization* quaternion_parameterization =
                        new QuaternionParameterization;
                    problem.SetParameterization(states[i].state, quaternion_parameterization);
                }
            }
            reportProgress();
            getchar();
            return true;
        }

#if 1

        virtual void reportProgress() {
            FILE *fx = fopen("X","w");
            FILE *fa = fopen("Ap","w"), *fb = fopen("Bp","w");
            fprintf(fx,"%e %e %e %e\n",common_parameters[0],common_parameters[1],
                    common_parameters[2],common_parameters[3]);
            for (size_t i=0;i<states.size();i++) {
                fprintf(fx,"%e %e %e %e %e\n",states[i].state[0],states[i].state[1],
                        states[i].state[2],states[i].state[3],states[i].depth);
                double As[3], Ap[3], Bs[3], Bp[3];
                double *Bscale = common_parameters;
                double k_depth = common_parameters[3];
                for (size_t j=0;j<3;j++) {
                    As[j] = states[i].a[j];
                    Bs[j] = states[i].m[j] * Bscale[j];
                }
                As[0] -= states[i].state[3]; // remove propulsion
                AngleAxisRotatePoint(states[i].state,As,Ap);
                Ap[2] -= k_depth*states[i].depth;
                AngleAxisRotatePoint(states[i].state,Bs,Bp);

                double t = (states[i].timestamp - states[0].timestamp)*24*3600;
                fprintf(fa,"%e %e %e %e\n",t,Ap[0],Ap[1],Ap[2]);
                fprintf(fb,"%e %e %e %e\n",t,Bp[0],Bp[1],Bp[2]);
            }
            fclose(fa); fclose(fb); fclose(fx);
            G.plot("set terminal x11 0");
            G.plot("plot \"Ap\" u 1:2 w l, \"Ap\" u 1:3 w l, \"Ap\" u 1:4 w l");
            G.plot("set terminal x11 1");
            G.plot("plot \"Bp\" u 1:2 w l, \"Bp\" u 1:3 w l, \"Bp\" u 1:4 w l");
        }
#endif
};


int main(int argc, char *argv[])
{
    google::ParseCommandLineFlags(&argc, &argv, true);
    google::InitGoogleLogging(argv[0]);
    if (FLAGS_input.empty()) {
        LOG(ERROR) << "Usage: bundle_adjustment_example --input=bal_problem";
        return 1;
    }



    OptimiseRotation problem;

    problem.load(FLAGS_input.c_str());

    problem.optimise();
    getchar();
}


