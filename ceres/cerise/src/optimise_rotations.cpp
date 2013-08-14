

#include <stdlib.h>
#include <stdio.h>
#include "cgnuplot/CGnuplot.h"
#include "ceres/ceres.h"
#include "gflags/gflags.h"
#include "glog/logging.h"
#include "rotation_errors.h"
#include "states.h"

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
    class OptimiseRotation {
        protected:
            cgnuplot::CGnuplot G;
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
            std::vector<DataLine> lines; 
            OptimisedOrientationSequence oos;


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
                lines.clear();
                FILE * fp = fopen(filename,"r");
                while (!feof(fp)) {
                    if ((lineLimit>0) && ((signed)lines.size() >= lineLimit)) {
                        break;
                    }
                    char line[4096] = {0,};
                    std::vector<double> dline;
                    if (fgets(line,4095, fp)!= NULL) {
                        DataLine dl;
                        if (dl.load(line)) {
                            lines.push_back(dl);
                        }
                    }
                }
                fclose(fp);
                printf("Loaded %d lines\n",(int)lines.size());

                LocalParameterization* quaternion_parameterization = NULL;
                oos.initialise(filename, lines);
                if (oos.use_quaternions && FLAGS_use_local_parameterization) {
                    quaternion_parameterization = new QuaternionParameterization;
                }
                for (size_t i=0;i<lines.size();i++) {
                    DataLine & dl(lines[i]);
                    OptimisedOrientation & oo(oos.states[i]);

                    LossFunction* loss_function;
                    CostFunction *cost_function;
                    // Acceleration
                    loss_function = FLAGS_robustify ? new HuberLoss(1.0) : NULL;
                    if (oos.use_quaternions) {
                        cost_function = new AutoDiffCostFunction<cerise::AccelerometerErrorQuat,3,4,4,1>(
                                new cerise::AccelerometerErrorQuat(dl.depth,dl.a[0],dl.a[1],dl.a[2], 10.0));
                    } else {
                        cost_function = new AutoDiffCostFunction<cerise::AccelerometerError,3,4,3,1>(
                                new cerise::AccelerometerError(dl.depth,dl.a[0],dl.a[1],dl.a[2], 10.0));
                    }
                    problem.AddResidualBlock(cost_function,loss_function,oos.common_parameters, 
                            oo.rotation,oo.propulsion);

                    // Magnetic field
                    loss_function = FLAGS_robustify ? new HuberLoss(100.0) : NULL;
                    if (oos.use_quaternions) {
                        cost_function = new AutoDiffCostFunction<cerise::MagnetometerErrorQuat,3,4,4>(
                                new cerise::MagnetometerErrorQuat(dl.m[0],dl.m[1],dl.m[2], 0.01));
                    } else {
                        cost_function = new AutoDiffCostFunction<cerise::MagnetometerError,3,4,3>(
                                new cerise::MagnetometerError(dl.m[0],dl.m[1],dl.m[2], 0.01));
                    }
                    problem.AddResidualBlock(cost_function,loss_function,oos.common_parameters, oo.rotation);

                    if (i>0) {
                        // Smoothness constraint
                        cost_function = new AutoDiffCostFunction<cerise::SmoothnessConstraint,1,1,1>(
                                new cerise::SmoothnessConstraint(1e2));
                        problem.AddResidualBlock(cost_function,NULL,oos.states[i-1].propulsion, oo.propulsion);
                    }
                    if (oos.use_quaternions && FLAGS_use_local_parameterization) {
                        problem.SetParameterization(oo.rotation, quaternion_parameterization);
                    }
                }
                reportProgress();
                getchar();
                return true;
            }

#if 1

            virtual void reportProgress() {
                oos.save("X");
                FILE *fa = fopen("Ap","w"), *fb = fopen("Bp","w");
                for (size_t i=0;i<lines.size();i++) {
                    DataLine & dl(lines[i]);
                    OptimisedOrientation & oo(oos.states[i]);
                    double As[3], Ap[3], Bs[3], Bp[3];
                    double t = (dl.timestamp - lines[0].timestamp)*24*3600;
                    if (oos.use_quaternions) {
                        for (size_t j=0;j<3;j++) {
                            As[j] = dl.a[j];
                            Bs[j] = dl.m[j] * oos.Bscale[j];
                        }
                        As[0] -= oo.propulsion[0]; // remove propulsion
                        QuaternionRotatePoint(oo.rotation,As,Ap);
                        Ap[2] -= oos.Kdepth[0]*dl.depth;
                        QuaternionRotatePoint(oo.rotation,Bs,Bp);
                    } else {
                        for (size_t j=0;j<3;j++) {
                            As[j] = dl.a[j];
                            Bs[j] = dl.m[j] * oos.Bscale[j];
                        }
                        As[0] -= oo.propulsion[0]; // remove propulsion
                        AngleAxisRotatePoint(oo.rotation,As,Ap);
                        Ap[2] -= oos.Kdepth[0]*dl.depth;
                        AngleAxisRotatePoint(oo.rotation,Bs,Bp);
                    }

                    fprintf(fa,"%e %e %e %e\n",t,Ap[0],Ap[1],Ap[2]);
                    fprintf(fb,"%e %e %e %e\n",t,Bp[0],Bp[1],Bp[2]);
                }
                fclose(fa); fclose(fb); 
                G.plot("set terminal wxt 0");
                G.plot("plot \"Ap\" u 1:2 w l, \"Ap\" u 1:3 w l, \"Ap\" u 1:4 w l");
                G.plot("set terminal wxt 1");
                G.plot("plot \"Bp\" u 1:2 w l, \"Bp\" u 1:3 w l, \"Bp\" u 1:4 w l");
                G.plot("set terminal wxt 3");
                if (oos.use_quaternions) {
                    G.plot("plot \"X\" u 0:6 w l, \"%s\" u 0:($9/500) w l",oos.input_file.c_str());
                } else {
                    G.plot("plot \"X\" u 0:5 w l, \"%s\" u 0:($9/500) w l",oos.input_file.c_str());
                }
            }
#endif
    };
};


int main(int argc, char *argv[])
{
    google::ParseCommandLineFlags(&argc, &argv, true);
    google::InitGoogleLogging(argv[0]);
    if (FLAGS_input.empty()) {
        LOG(ERROR) << "Usage: optimise_rotation --input=preload.txt" << endl << "\twhere preload has been prepared by the prepare.m matlab script\n";
        return 1;
    }



    cerise::OptimiseRotation problem;

    problem.load(FLAGS_input.c_str());

    problem.optimise();
    getchar();
}


