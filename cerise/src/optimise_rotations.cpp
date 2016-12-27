

#include <stdlib.h>
#include <stdio.h>
#include "cgnuplot/CGnuplot.h"
#include "ceres/ceres.h"
#include "gflags/gflags.h"
#include "glog/logging.h"
#include "rotation_errors.h"
#include "states.h"

DEFINE_string(input, "", "Input File name");
DEFINE_string(output, "", "Input File name");
DEFINE_string(gps, "", "GPS File name");
DEFINE_bool(raw, false, "Use a raw input file (11 columns)");
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
DEFINE_bool(interactive, false, "Wait for user key presses.");
DEFINE_bool(display, false, "Display plot during progress");

DEFINE_double(ftol, 1e-7, "Function tolerance");
DEFINE_double(eta, 1e-2, "Default value for eta. Eta determines the "
             "accuracy of each linear solve of the truncated newton step. "
             "Changing this parameter can affect solve performance.");

DEFINE_bool(use_local_parameterization, false, "For quaternions, use a local "
            "parameterization.");

DEFINE_int32(num_lines, -1, "Number of data lines to consider.");
DEFINE_int32(line_step, -1, "Number of data lines to consider per iteration.");
DEFINE_int32(num_threads, 1, "Number of threads.");
DEFINE_int32(num_iterations, 100, "Number of iterations.");
DEFINE_double(max_solver_time, 1e32, "Maximum solve time in seconds.");
DEFINE_bool(nonmonotonic_steps, false, "Trust region algorithm can use"
            " nonmonotic steps.");
DEFINE_bool(mag,false,"Estimate the magnetic field vector");

DEFINE_string(solver_log, "", "File to record the solver execution to.");



using namespace ceres;

namespace cerise{ 
    class OptimiseRotation {
        protected:
            cgnuplot::CGnuplot G;
            class DisplayCallback: public ceres::IterationCallback { 
                protected:
                    OptimiseRotation & problem;
                    size_t startLine, endLine;
                public: 
                    DisplayCallback(OptimiseRotation & p) : problem(p) {}
                    void setLines(size_t start,size_t end) {
                        startLine = start;
                        endLine = end;
                    }
                    virtual ceres::CallbackReturnType operator()(const 
                            ceres::IterationSummary& summary) { 
                        problem.reportProgress(startLine,endLine,true);
                        // if (FLAGS_interactive) {
                        //     getchar();
                        // }
                        return ceres::SOLVER_CONTINUE;
                    } 
            };

            DisplayCallback display;
            std::vector<DataLine> lines; 
            typedef std::map<double,GPSLine> GPSMap; 
            GPSMap gpslines; 
            OptimisedOrientationSequence oos;
            std::string filename;


            Problem *problem;
            void SetLinearSolver(Solver::Options* options) {
                CHECK(StringToLinearSolverType(FLAGS_linear_solver,
                            &options->linear_solver_type));
                CHECK(StringToPreconditionerType(FLAGS_preconditioner,
                            &options->preconditioner_type));
                // CHECK(StringToSparseLinearAlgebraLibraryType(
                //             FLAGS_sparse_linear_algebra_library,
                //             &options->sparse_linear_algebra_library));
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
                options->use_inner_iterations = FLAGS_inner_iterations;
            }

            const GPSLine & getClosestGPSFix(double ts) {
                assert(gpslines.size()>0);
                GPSMap::const_iterator it2 = gpslines.upper_bound(ts);

                if (it2 == gpslines.begin()) {
                    return it2->second;
                }

                GPSMap::const_iterator it1 = it2;
                it1 --;
                if (it2 == gpslines.end())  {
                    return it1->second;
                }
                if (fabs(it1->first-ts) > fabs(it2->first-ts)) {
                    return it2->second;
                } else {
                    return it1->second;
                }
            }

        public:
            OptimiseRotation() : display(*this), problem(NULL) {}

            void optimise(size_t startLine, size_t endLine) {
                assert(problem);
                Solver::Options options;
                SetMinimizerOptions(&options);
                SetLinearSolver(&options);
                Solver::Summary summary;
                Solve(options, problem, &summary);
                std::cout << summary.FullReport() << "\n";
                if (!FLAGS_output.empty()) {
                    oos.save(FLAGS_output);
                }
                reportProgress(startLine,endLine);
                // sleep(3);
            }
            bool loadGPS(const char *gps, int lineLimit=-1) {
                gpslines.clear();
                FILE * fp = fopen(gps,"r");
                if (!fp) {
                    LOG(ERROR) << "Couldn't load gps file '" << gps << "'\n";
                    return false;
                }
                unsigned int target = 500000;
                while (!feof(fp)) {
                    char line[4096] = {0,};
                    if (fgets(line,4095, fp)!= NULL) {
                        if ((lineLimit>0) && ((signed)gpslines.size() >= lineLimit)) {
                            break;
                        }
                        GPSLine dl;
                        if (dl.load(line)) {
                            gpslines.insert(GPSMap::value_type(dl.timestamp,dl));
                        }
                        if (gpslines.size() >= target) {
                            printf("GPS: Loaded %d lines\n",(int)gpslines.size());
                            target += 500000;
                        }
                    }
                }
                fclose(fp);
                printf("Finished loading %d gps lines\n",(int)gpslines.size());
                assert(gpslines.size() >= 2);
                return true;
            }

            bool load(const char * _filename, const char *gps, int lineLimit=-1) {
                filename = _filename;
                lines.clear();
                FILE * fp = fopen(_filename,"r");
                unsigned int target = 500000;
                while (!feof(fp)) {
                    if ((lineLimit>0) && ((signed)lines.size() >= lineLimit)) {
                        break;
                    }
                    if (lines.size() >= target) {
                        printf("Data: Loaded %d lines\n",(int)lines.size());
                        target += 500000;
                    }
                        
                    char line[4096] = {0,};
                    std::vector<double> dline;
                    if (fgets(line,4095, fp)!= NULL) {
                        DataLine dl;
                        if (FLAGS_raw) {
                            if (dl.load_raw(line)) {
                                lines.push_back(dl);
                            }
                        } else {
                            if (dl.load(line)) {
                                lines.push_back(dl);
                            }
                        }
                    }
                }
                fclose(fp);
                printf("Finished loading %d data lines\n",(int)lines.size());
                // prepare a map between dive_ids and variable numbers
                std::map<size_t,size_t> dive_ids;
                for (size_t i=0;i<lines.size();i++) {
                    if (dive_ids.find(lines[i].dive_id) == dive_ids.end()) {
                        dive_ids[lines[i].dive_id] = dive_ids.size();
                    }
                }

#if 1
                if (!loadGPS(gps,lineLimit)) {
                    return false;
                }
#endif
                return true;
            }

            size_t numLines() {
                return lines.size();
            }

            bool createProblem(size_t beginLine, size_t endLine) {
                delete problem;
                problem = new Problem();
                display.setLines(beginLine,endLine);

                LocalParameterization* quaternion_parameterization = NULL;
                oos.initialise(filename, lines, beginLine, endLine);
                if (oos.use_quaternions && FLAGS_use_local_parameterization) {
                    quaternion_parameterization = new QuaternionParameterization;
                }
                for (size_t i=beginLine;i<std::min<size_t>(endLine,lines.size());i++) {
                    DataLine & dl(lines[i]);
                    const GPSLine & gps = getClosestGPSFix(dl.timestamp);
                    OptimisedOrientation & oo(oos.states[i-beginLine]);

                    LossFunction* loss_function;
                    CostFunction *cost_function;
                    // Acceleration
                    loss_function = FLAGS_robustify ? new HuberLoss(1.0) : NULL;
                    if (oos.use_quaternions) {
                        if (FLAGS_mag) {
                            cost_function = new AutoDiffCostFunction<cerise::AccelerometerErrorQuat,3,7,4,1>(
                                    new cerise::AccelerometerErrorQuat(dl.depth,dl.a[0],dl.a[1],dl.a[2], 100.0));
                        } else {
                            cost_function = new AutoDiffCostFunction<cerise::AccelerometerErrorQuat,3,4,4,1>(
                                    new cerise::AccelerometerErrorQuat(dl.depth,dl.a[0],dl.a[1],dl.a[2], 100.0));
                        }
                    } else {
                        cost_function = new AutoDiffCostFunction<cerise::AccelerometerError,3,4,3,1>(
                                new cerise::AccelerometerError(dl.depth,dl.a[0],dl.a[1],dl.a[2], 100.0));
                    }
                    problem->AddResidualBlock(cost_function,loss_function,oos.common_parameters, 
                            oo.rotation,oo.propulsion);

                    // Magnetic field
                    loss_function = FLAGS_robustify ? new HuberLoss(100.0) : NULL;
                    if (oos.use_quaternions) {
                        if (FLAGS_mag) {
                            cost_function = new AutoDiffCostFunction<cerise::MagnetometerFieldErrorQuat,3,7,4>(
                                    new cerise::MagnetometerFieldErrorQuat(dl.m[0],dl.m[1],dl.m[2], 0.001));
                        } else {
                            cost_function = new AutoDiffCostFunction<cerise::MagnetometerErrorQuat,3,4,4>(
                                    new cerise::MagnetometerErrorQuat(dl.m[0],dl.m[1],dl.m[2], 
                                        gps.B[0],gps.B[1],gps.B[2], 0.001));
                        }
                    } else {
                        cost_function = new AutoDiffCostFunction<cerise::MagnetometerError,3,4,3>(
                                new cerise::MagnetometerError(dl.m[0],dl.m[1],dl.m[2], 
                                    gps.B[0],gps.B[1],gps.B[2], 0.001));
                    }
                    problem->AddResidualBlock(cost_function,loss_function,oos.common_parameters, oo.rotation);

                    if (i-beginLine>0) {
                        // Smoothness constraint
                        cost_function = new AutoDiffCostFunction<cerise::SmoothnessConstraint,1,1,1>(
                                new cerise::SmoothnessConstraint(5e1));
                        problem->AddResidualBlock(cost_function,NULL,oos.states[i-beginLine-1].propulsion, oo.propulsion);
                    }
                    if (oos.use_quaternions && FLAGS_use_local_parameterization) {
                        problem->SetParameterization(oo.rotation, quaternion_parameterization);
                    }
                }
                if (FLAGS_interactive) {
                    reportProgress(beginLine,endLine);
                    printf("Press enter to start\n");
                    getchar();
                }
                return true;
            }

            void saveState(long long unsigned int startLine, long long unsigned int endLine) {
                char buffer[1024];
                sprintf(buffer,"X%08llu",startLine);
                oos.save(buffer);
                sprintf(buffer,"Ap%08llu",startLine);
                FILE *fa = fopen(buffer,"w");
                sprintf(buffer,"Bp%08llu",startLine);
                FILE *fb = fopen(buffer,"w");
                for (size_t i=startLine;i<endLine;i++) {
                    DataLine & dl(lines[i]);
                    OptimisedOrientation & oo(oos.states[i-startLine]);
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
            }

#if 1

            virtual void reportProgress(long long unsigned int startLine, long long unsigned int endLine, bool intermediate=false) {
                oos.save("X");
                FILE *fa = fopen("Ap","w"), *fb = fopen("Bp","w");
                for (size_t i=startLine;i<endLine;i++) {
                    DataLine & dl(lines[i]);
                    OptimisedOrientation & oo(oos.states[i-startLine]);
                    double As[3], Ap[3], Bs[3], Bp[3];
                    double t = (dl.timestamp - lines[startLine].timestamp)*24*3600;
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
                if (FLAGS_display) {
                    G.plot("set terminal x11 1;set grid");
                    G.plot("plot \"Bp\" u 1:2 w l, \"Bp\" u 1:3 w l, \"Bp\" u 1:4 w l");
                    if (!intermediate) {
                        G.plot("set terminal x11 3;set grid");
                        int n = oos.states.size();
                        int depth_col=9;
                        if (FLAGS_raw) {
                            depth_col=3;
                        }
                        if (oos.use_quaternions) {
                            G.plot("plot [0:%d] \"X\" u 0:6 w l, \"%s\" u 0:($%d/500) w l",n,oos.input_file.c_str(),depth_col);
                        } else {
                            G.plot("plot [0:%d] \"X\" u 0:5 w l, \"%s\" u 0:($%d/500) w l",n,oos.input_file.c_str(),depth_col);
                        }
                    } 
                    G.plot("set terminal x11 0;set grid");
                    G.plot("plot \"Ap\" u 1:2 w l, \"Ap\" u 1:3 w l, \"Ap\" u 1:4 w l");
                }
            }
#endif
    };
};


int main(int argc, char *argv[])
{
    gflags::ParseCommandLineFlags(&argc, &argv, true);
    google::InitGoogleLogging(argv[0]);
    if (FLAGS_input.empty()) {
        LOG(ERROR) << "Usage: optimise_rotation --input=preload.txt\n\twhere preload has been prepared by the prepare.m matlab script\n";
        return 1;
    }



    cerise::OptimiseRotation problem;

    problem.load(FLAGS_input.c_str(),FLAGS_gps.c_str(),FLAGS_num_lines);
    if (FLAGS_line_step < 0) {FLAGS_line_step = problem.numLines();}
    if (FLAGS_num_lines < 0) {FLAGS_num_lines = problem.numLines();}
    FLAGS_num_lines = std::min<size_t>(problem.numLines(),FLAGS_num_lines);

    for (size_t startLine=0;startLine < (size_t)FLAGS_num_lines; startLine += FLAGS_line_step) { 
        size_t lastLine = std::min<size_t>(FLAGS_num_lines,startLine+FLAGS_line_step);
        printf("Starting optimisation in [%d,%d]\n",(int)startLine,(int)lastLine);
        problem.createProblem(startLine,lastLine);
        problem.optimise(startLine,lastLine);
        problem.saveState(startLine,lastLine);
        if (FLAGS_interactive) {
            printf("Optimisation [%d,%d] complete\n",(int)startLine,(int)lastLine);
            getchar();
        }
    }

    return 0;
}


