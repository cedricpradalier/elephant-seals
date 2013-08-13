#ifndef OPTIMISATION_PROBLEM_H
#define OPTIMISATION_PROBLEM_H

#include <boost/shared_ptr.hpp>
#include "elephants/ErrorTerm.h"

template <class Common, class State, class Observation>
class OptimisationProblem
{
    protected:
        typedef ObservationError<Common,State,Observation> ObservationError_;
        typedef UnassociatedObservationError<Common,State,Observation> UnassociatedObservationError_;
        typedef TransitionError<Common,State,Observation> TransitionError_;
        typedef UnassociatedTransitionError<Common,State,Observation> UnassociatedTransitionError_;

        std::vector< boost::shared_ptr<ObservationError_> > observationErrors;
        std::vector< boost::shared_ptr<UnassociatedObservationError_> > unassObservationErrors;
        std::vector< boost::shared_ptr<TransitionError_> > transitionErrors;
        std::vector< boost::shared_ptr<UnassociatedTransitionError_> > unassTransitionErrors;

        size_t stateSize, nEq;
        Common common, candidateCommon;
        std::vector< State > states, candidateStates;
        std::vector< Observation > observations;
        std::vector<size_t> stateIndex; 


        typedef Eigen::SparseMatrix<float> SpMat; // declares a column-major sparse matrix type of double

        Eigen::VectorXf F;
        SpMat W, J;

    public:
        OptimisationProblem() {
        }

        void updateMatrices() {
            stateSize = common.getDimension();
            for (size_t i=0;i<states.size();i++) {
                stateIndex[i] = stateSize;
                stateSize += states[i].getDimension();
            }
            nEq = 0;
            for (size_t i=0;i<observationErrors.size();i++) {
                nEq += observationErrors[i]->getDimension();
            }
            for (size_t i=0;i<unassObservationErrors.size();i++) {
                nEq += unassObservationErrors[i]->getDimension();
            }
            for (size_t i=0;i<transitionErrors.size();i++) {
                nEq += transitionErrors[i]->getDimension();
            }
            for (size_t i=0;i<transitionErrors.size();i++) {
                nEq += unassTransitionErrors[i]->getDimension();
            }
            F = Eigen::VectorXf(nEq);
            W = SpMat(nEq,nEq);
            J = SpMat(nEq,stateSize);
        }

        void computeCandidateState(const Eigen::VectorXf & delta) {
            candidateCommon = common;
            candidateStates = states;
            candidateCommon.step(delta.block(0,common.getDimension()));
            for (size_t i=0;i<states.size();i++) {
                candidateStates[i].step(delta.block(stateIndex[i],states[i].getDimension()));
            }
        }

        void optimise() {
            double mu = 1.0;
            double scale = 1.1;
            SpMat I(stateSize,stateSize);
            I.setIdentity();
            while (mu < 20.0) {
                updateError();
                SpMat JW = J.transpose() * W;
                double currentNorm = F.norm();

                while (mu < 20.0) {
                    // Solving:
                    SpMat A = JW * J + mu * I;
                    Eigen::SimplicialCholesky<SpMat> chol(A);  // performs a Cholesky factorization of A
                    Eigen::VectorXf delta = chol.solve(JW * F);         // use the factorization to solve for the given right hand side
                    computeCandidateState(delta);
                    double newNorm = computeCandidateNorm();
                    if (newNorm < currentNorm) {
                        states = candidateStates;
                        common = candidateCommon;
                        mu /= scale;
                        break;
                    } else {
                        mu *= scale;
                    }
                }
            }

        }


        void computeCandidateNorm() {
            size_t eqIndex = 0;
            Eigen::VectorXf F(nEq);
            for (size_t k=0;k<observationErrors.size();k++) {
                shared_ptr<ObservationError_> oe = observationErrors[k];
                for (size_t i=0;i<states.size();i++) {
                    oe->update(common,states[i],observations[i]);
                    size_t neq = oe->getDimension();
                    if (neq) {
                        F.block(eqIndex,neq) = oe->getError();
                    }
                }
                eqIndex += oe->getDimension();
            }
            for (size_t k=0;k<unassObservationErrors.size();k++) {
                shared_ptr<UnassociatedObservationError_> oe = unassObservationErrors[k];
                int res = oe->update(common,states);
                size_t neq = oe->getDimension();
                if (neq) {
                    F.block(eqIndex,neq) = oe->getError();
                }
                eqIndex += oe->getDimension();
            }
            for (size_t k=0;k<transitionErrors.size();k++) {
                shared_ptr<TransitionError_> te = transitionErrors[k];
                for (size_t i=0;i<states.size()-1;i++) {
                    te->update(common,states[i],observations[i],states[i+1],observations[i+1]);
                    size_t neq = te->getDimension();
                    if (neq) {
                        F.block(eqIndex,neq) = te->getError();
                    }
                }
                eqIndex += oe->getDimension();
            }
            for (size_t k=0;k<unassTransitionErrors.size();k++) {
                shared_ptr<UnassociatedTransitionError_> te = unassTransitionErrors[k];
                std::pair<int,int> res = te->update(common,states);
                size_t neq = te->getDimension();
                if (neq) {
                    F.block(eqIndex,neq) = oe->getError();
                }
                eqIndex += oe->getDimension();
            }

            return F.norm();
        }

        void updateError() {
            size_t ncom = common.getDimension();
            size_t eqIndex = 0;
            F = Eigen::VectorXf(nEq);
            if (!valueOnly) {
                W = SpMat(nEq,nEq);
                J = SpMat(nEq,stateSize);
            }
            for (size_t k=0;k<observationErrors.size();k++) {
                shared_ptr<ObservationError_> oe = observationErrors[k];
                for (size_t i=0;i<states.size();i++) {
                    oe->update(common,states[i],observations[i]);
                    size_t neq = oe->getDimension();
                    if (neq) {
                        F.block(eqIndex,neq) = oe->getError();
                        W.block(eqIndex,eqIndex,neq,neq) = oe->getWeights();
                        J.block(eqIndex,stateIndex[i],neq,states[i].getDimension()) = oe->getJacobian();
                        if (ncom) {
                            J.block(eqIndex,0,neq,ncom) = oe->getJacobianWrtCommonState();
                        }
                    }
                }
                eqIndex += oe->getDimension();
            }
            for (size_t k=0;k<unassObservationErrors.size();k++) {
                shared_ptr<UnassociatedObservationError_> oe = unassObservationErrors[k];
                int res = oe->update(common,states);
                size_t neq = oe->getDimension();
                if (neq) {
                    F.block(eqIndex,neq) = oe->getError();
                    W.block(eqIndex,eqIndex,neq,neq) = oe->getWeights();
                    J.block(eqIndex,stateIndex[res],neq,states[res].getDimension()) = oe->getJacobian();
                    if (ncom) {
                        J.block(eqIndex,0,neq,ncom) = oe->getJacobianWrtCommonState();
                    }
                }
                eqIndex += oe->getDimension();
            }
            for (size_t k=0;k<transitionErrors.size();k++) {
                shared_ptr<TransitionError_> te = transitionErrors[k];
                for (size_t i=0;i<states.size()-1;i++) {
                    te->update(common,states[i],observations[i],states[i+1],observations[i+1]);
                    size_t neq = te->getDimension();
                    if (neq) {
                        F.block(eqIndex,neq) = te->getError();
                        W.block(eqIndex,eqIndex,neq,neq) = te->getWeights();
                        J.block(eqIndex,stateIndex[i],neq,
                                states[i].getDimension()+states[i+1].getDimension()) = te->getJacobian();
                        if (ncom) {
                            J.block(eqIndex,0,neq,ncom) = te->getJacobianWrtCommonState();
                        }
                    }
                }
                eqIndex += oe->getDimension();
            }
            for (size_t k=0;k<unassTransitionErrors.size();k++) {
                shared_ptr<UnassociatedTransitionError_> te = unassTransitionErrors[k];
                std::pair<int,int> res = te->update(common,states);
                size_t neq = te->getDimension();
                if (neq) {
                    F.block(eqIndex,neq) = oe->getError();
                    W.block(eqIndex,eqIndex,neq,neq) = oe->getWeights();
                    J.block(eqIndex,stateIndex[res.first],neq,states[res.first].getDimension()) = 
                        te->getJacobian().block(0,0,neq,states[res.first].getDimension());
                    J.block(eqIndex,stateIndex[res.second],neq,states[res.second].getDimension()) = 
                        te->getJacobian().block(0,states[res.first].getDimension(),neq,states[res.second].getDimension());
                    if (ncom) {
                        J.block(eqIndex,0,neq,ncom) = oe->getJacobianWrtCommonState();
                    }
                }
                eqIndex += oe->getDimension();
            }
        }

};

#endif // OPTIMISATION_PROBLEM_H
