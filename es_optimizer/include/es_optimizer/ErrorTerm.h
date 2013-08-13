#ifndef ERROR_TERM_H
#define ERROR_TERM_H

#include <Eigen/Core.h>


class StateDescription {
    public:
        StateDescription() {}
        ~virtual StateDescription() {}
        virtual size_t getDimension() = 0;

        virtual void step(const Eigen::VectorXf & delta);
};

class CommonState : public StateDescription{
    public:
        CommonState() {}
        virtual ~CommonState() {}
};


class ObservationDescription {
    public:
        ObservationDescription() {}
};

class ErrorTerm {

    public:
        ErrorTerm() {}

        virtual ~ErrorTerm() {}

        virtual size_t getDimension() = 0;
        virtual const Eigen::VectorXf & getError() = 0;
        virtual const Eigen::MatrixXf & getWeights() = 0;
        virtual const Eigen::MatrixXf & getJacobian() = 0;
        virtual const Eigen::MatrixXf & getJacobianWrtCommonState() = 0;

};

template <class Common, class State, class Observation>
class ObservationError : public ErrorTerm {
    public:
        ObservationError() : ErrorTerm() {}
        virtual ~ObservationError() {}

        virtual void update(const Common & c, const State & s, const Observation & o) = 0;
};

template <class Common, class State, class Observation>
class UnassociatedObservationError : public ErrorTerm {
    public:
        UnassociatedObservationError() : ErrorTerm() {}

        // Returns the index of the associated state
        virtual size_t update(const Common & c, const std::vector<State> & s) = 0;
};


template <class Common, class State, class Observation>
class TransitionError : public ErrorTerm {
    public:
        TransitionError() : ErrorTerm() {}

        virtual void update(const Common & c, const State & s1, const Observation & o1, const State & s2, const Observation & o2) = 0;
};

template <class Common, class State, class Observation>
class UnassociatedTransitionError : public ErrorTerm {
    public:
        UnassociatedTransitionError() : ErrorTerm() {}

        virtual std::pair<size_t,size_t> update(const Common & c, const std::vector<State> & state ) = 0;
};




#endif // ERROR_TERM_H
