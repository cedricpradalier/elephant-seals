
#include <stdlib.h>
#include <stdio.h>
#include <ros/ros.h>
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
DEFINE_bool(interactive, false, "Wait for user key presses.");
DEFINE_bool(display, false, "Display plot during progress");

DEFINE_double(ftol, 1e-7, "Function tolerance");
DEFINE_double(eta, 1e-2, "Default value for eta. Eta determines the "
             "accuracy of each linear solve of the truncated newton step. "
             "Changing this parameter can affect solve performance.");

DEFINE_int32(num_lines, -1, "Number of data lines to consider.");
DEFINE_int32(num_threads, 1, "Number of threads.");
DEFINE_int32(num_iterations, 100, "Number of iterations.");
DEFINE_double(max_solver_time, 1e32, "Maximum solve time in seconds.");
DEFINE_bool(nonmonotonic_steps, false, "Trust region algorithm can use"
            " nonmonotic steps.");
DEFINE_bool(mag,false,"Estimate the magnetic field vector");

DEFINE_string(solver_log, "", "File to record the solver execution to.");



using namespace ceres;

namespace cerise{ 
    struct DataLine {
        double timestamp;
        double X[3]; // ground truth position
        double V[3]; // ground truth velocity
        double rpy[3]; // Euler angles in rad
        // Measurements
        double rpy_init[3]; // Euler angles in rad
        double a[3];
        double m[3];

        bool load(const std::string & line) {
            char buffer[line.size()];
            double f=0.0;
            int n=0;
            const char * lp = line.c_str();
            std::vector<double> dline;
            while (sscanf(lp, " %s%n", buffer, &n)==1) {
                if (sscanf(buffer,"%le",&f)==1) {
                    dline.push_back(f);
                } else {
                    dline.push_back(NAN);
                }
                lp += n;
            }
            if (dline.size() < 16) {
                LOG(ERROR) << "Not enough data filed on input line\n";
                return false;
            }
            timestamp = dline[0];
            for (int i=0;i<3;i++) { 
                X[i] = dline[1+i]; 
                V[i] = dline[4+i]; 
                rpy[i] = dline[7+i];
                a[i] = dline[10+i];
                m[i] = dline[13+i];
            }
            double na = ::sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
            double nm = ::sqrt(m[0]*m[0]+m[1]*m[1]+m[2]*m[2]);

            for (int i=0;i<3;i++) { 
                a[i] /= na;
                m[i] /= nm;
            }
            assert(fabs(na-10.0) < 2.0);
            rpy_init[1] = ::asin(a[0] / na);
            rpy_init[0] = -::asin((a[1]/na)/::cos(rpy_init[1]));
            rpy_init[2] = 0;

            return true;
        }
    };

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
                        problem.reportProgress(true);
                        // if (FLAGS_interactive) {
                        //     getchar();
                        // }
                        return ceres::SOLVER_CONTINUE;
                    } 
            };

            DisplayCallback display;
            std::vector<DataLine> lines; 
            boost::shared_ptr<double> P;
            boost::shared_ptr<double> Q;


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
                options->function_tolerance = FLAGS_ftol;
                options->max_solver_time_in_seconds = FLAGS_max_solver_time;
                options->use_nonmonotonic_steps = FLAGS_nonmonotonic_steps;
                if (FLAGS_display) {
                    options->callbacks.push_back(&display);
                }
                options->update_state_every_iteration = true; 

                CHECK(StringToTrustRegionStrategyType(FLAGS_trust_region_strategy,
                            &options->trust_region_strategy_type));
                CHECK(StringToDoglegType(FLAGS_dogleg, &options->dogleg_type));
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
                std::cout << "P: " << P.get()[0] 
                    << " " << P.get()[1]
                    << " " << P.get()[2];
                reportProgress();
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
                        if (line[0]=='#') continue;
                        DataLine dl;
                        if (dl.load(line)) {
                            lines.push_back(dl);
                        }
                    }
                }
                fclose(fp);
                printf("Loaded %d lines\n",(int)lines.size());

                LocalParameterization* quaternion_parameterization = NULL;
                P.reset(new double[2]);
                P.get()[0] = 1.0;
                P.get()[1] = 1.0;
                P.get()[2] = 1.0;
                Q.reset(new double[4*lines.size()]);
                quaternion_parameterization = new QuaternionParameterization;
                for (size_t i=0;i<lines.size();i++) {
                    DataLine & dl(lines[i]);
                    double *q = Q.get()+4*i;
                    double mat[9] = {1,0,0,0,1,0,0,0,1};
                    double log[3];
                    EulerAnglesToRotationMatrix<double>(dl.rpy_init,3,mat);
                    RotationMatrixToAngleAxis<double>(mat,log);
                    AngleAxisToQuaternion<double>(log,q);
                    LossFunction* loss_function;
                    CostFunction *cost_function;
                    // Acceleration
                    loss_function = FLAGS_robustify ? new HuberLoss(1.0) : NULL;
                    cost_function = new AutoDiffCostFunction<cerise::AccelerometerErrorQuat,3,3,4>(
                            new cerise::AccelerometerErrorQuat(dl.a[0],dl.a[1],dl.a[2], 1.0));
                    problem.AddResidualBlock(cost_function,loss_function,P.get(), Q.get()+4*i);

                    // Magnetic field
                    loss_function = FLAGS_robustify ? new HuberLoss(100.0) : NULL;
                    cost_function = new AutoDiffCostFunction<cerise::MagnetometerErrorQuat,3,3,4>(
                            new cerise::MagnetometerErrorQuat(dl.m[0],dl.m[1],dl.m[2], 1.0));
                    problem.AddResidualBlock(cost_function,loss_function,P.get(), Q.get()+4*i);

                    problem.SetParameterization(Q.get()+4*i, quaternion_parameterization);
                }
                if (FLAGS_interactive) {
                    reportProgress();
                    getchar();
                }
                return true;
            }

#if 1
            virtual void reportProgress(bool intermediate=false) {
                FILE *fa = fopen("Ap","w"), *fb = fopen("Bp","w");
                for (size_t i=0;i<lines.size();i++) {
                    DataLine & dl(lines[i]);
                    double As[3], Ap[3], Bs[3], Bp[3];
                    double t = dl.timestamp;
                    for (size_t j=0;j<3;j++) {
                        As[j] = dl.a[j] * P.get()[0];
                        Bs[j] = dl.m[j] * P.get()[1];
                    }
                    QuaternionRotatePoint(Q.get()+4*i,As,Ap);
                    QuaternionRotatePoint(Q.get()+4*i,Bs,Bp);

                    fprintf(fa,"%e %e %e %e\n",t,Ap[0],Ap[1],Ap[2]);
                    fprintf(fb,"%e %e %e %e\n",t,Bp[0],Bp[1],Bp[2]);
                }
                fclose(fa); fclose(fb); 
                if (FLAGS_display) {
                    G.plot("set terminal x11 1;set grid");
                    G.plot("plot \"Bp\" u 1:2 w l, \"Bp\" u 1:3 w l, \"Bp\" u 1:4 w l");
                    G.plot("set terminal x11 0;set grid");
                    G.plot("plot \"Ap\" u 1:2 w l, \"Ap\" u 1:3 w l, \"Ap\" u 1:4 w l");
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

    problem.load(FLAGS_input.c_str(),FLAGS_num_lines);

    problem.optimise();
    if (FLAGS_interactive) {
        getchar();
    }

    return 0;
}


