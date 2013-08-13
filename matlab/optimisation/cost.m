function [F,W,J] = cost(common,state,obs,problem)
% common is a struct containing the values of the common variables
% state is an array of state struct
% obs is an array of observation struct
% problem is a struct with the following 5 members:
%   state_dof is the number of degrees of freedom of the state we're optimising
%   The following 4 members are cell arrays of functions
%   [f,w,jc,js] = f_state contains function of a single state f(common, s,o) where s is a state
%   struct and o an observation struct
%   [f,w,jc,js] = f_delta contains function of two subsequent states f(common, s1,s2,o1,o2) where s1,s2 are state
%   structs and o1,o2 are observation structs
%   [f,w,jc,js,s] = f_async is a function of a single state f(common, state) used for
%   asynchronous events, s in the output is the index of the state being
%   measured
%   [f,w,jc,j1,s1,j2,s2] = f_async is a function of a two states f(common, state) used for
%   asynchronous events linking two states. s1 and s2 in the output are the index of the state being
%   measured, and j1, j2 the respective jacobians.
% In all the functions, f is the measure of the cost to be optimised, w is the
% weight associated with this particular measurement, jc (resp js) is the jacobian of f
% with respect to the common state (resp the state or states). If a function thinks it is not relevant
% at a given state, it can just returns f = []

    N = size(state,2);
    F = [];
    J = [];
    W = [];
    if nargout > 1
        W=sparse(0,0);
    end
    if nargout > 2
        J=sparse(0,0);
    end
    eq = 1;
    disp '   f_state'
    for k=1:size(problem.f_state,2)
        sprintf('      %d\n',k)
        for i=1:N
            vbase = problem.common_dof + (i-1)*problem.state_dof;
            if nargout > 2
                [f,w,jc,js] = problem.f_state{k}(common,state{i},obs{i});
            elseif nargout > 1
                [f,w] = problem.f_state{k}(common,state{i},obs{i});
            else
                f = problem.f_state{k}(common,state{i},obs{i});
            end
            n_eq = size(f,1);
            if n_eq==0
                continue
            end
            F(eq:eq+n_eq-1) = f;
            if nargout > 1
                assert(isreal(jc));
                assert(isreal(js));
                if problem.common_dof
                    J(eq:eq+n_eq-1, 1:problem.common_dof) = jc;
                end
                J(eq:eq+n_eq-1, vbase+1:vbase+problem.state_dof) = js;
            end
            if nargout > 2
                W(eq:eq+n_eq-1,eq:eq+n_eq-1) = w;
            end
            eq = eq+n_eq;
        end
    end

    disp '   f_delta'
    for k=1:size(problem.f_delta,2)
        sprintf('      %d\n',k)
        for i=2:N
            vbaseprev = problem.common_dof + (i-2)*problem.state_dof;
            if nargout > 2
                [f,w,jc,js] = problem.f_delta{k}(common,state{i-1},state{i},obs{i-1},obs{i});
            elseif nargout > 1
                [f,w] = problem.f_delta{k}(common,state{i-1},state{i},obs{i-1},obs{i});
            else
                f = problem.f_delta{k}(common,state{i-1},state{i},obs{i-1},obs{i});
            end
            n_eq = size(f,1);
            if n_eq==0
                continue
            end
            F(eq:eq+n_eq-1) = f;
            if nargout > 1
                assert(isreal(jc));
                assert(isreal(js));
                if problem.common_dof
                    J(eq:eq+n_eq-1, 1:problem.common_dof) = jc;
                end
                J(eq:eq+n_eq-1, vbaseprev+1:vbaseprev+2*problem.state_dof) = js;
            end
            if nargout > 2
                W(eq:eq+n_eq-1,eq:eq+n_eq-1) = w;
            end
            eq = eq+n_eq;
        end
    end

    disp '   f_async'
    for k=1:size(problem.f_async,2)
        sprintf('      %d\n',k)
        if nargout > 2
            [f,w,jc,j,s] = problem.f_async{k}(common,state);
            vbase = problem.common_dof + (s-1)*problem.state_dof;
        elseif nargout > 1
            [f,w] = problem.f_async{k}(common,state);
        else 
            f = problem.f_async{k}(common,state);
        end
        n_eq = size(f,1);
        if n_eq==0
            continue
        end
        F(eq:eq+n_eq-1) = f;
        if nargout > 1
            assert(isreal(jc));
            assert(isreal(j));
            if problem.common_dof
                J(eq:eq+n_eq-1, 1:problem.common_dof) = jc;
            end
            J(eq:eq+n_eq-1, vbase+1:vbase+problem.state_dof) = j;
        end
        if nargout > 2
            W(eq:eq+n_eq-1,eq:eq+n_eq-1) = w;
        end
        eq = eq+n_eq;
    end

    disp '   f_async_2'
    for k=1:size(problem.f_async2,2)
        sprintf('      %d\n',k)
        if nargout > 2
            [f,w,jc,j1,s1,j2,s2] = problem.f_async{k}(common,state);
            vbase1 = problem.common_dof + (s1-1)*problem.state_dof;
            vbase2 = problem.common_dof + (s2-1)*problem.state_dof;
        elseif nargout > 1
            [f,w] = problem.f_async{k}(common,state);
        else
            f = problem.f_async{k}(common,state);
        end
        n_eq = size(f,1);
        if n_eq==0
            continue
        end
        F(eq:eq+n_eq-1) = f;
        if nargout > 1
            assert(isreal(jc));
            assert(isreal(j1));
            assert(isreal(j2));
            if problem.common_dof
                J(eq:eq+n_eq-1, 1:problem.common_dof) = jc;
            end
            J(eq:eq+n_eq-1, vbase1+1:vbase1+problem.state_dof) = j1;
            J(eq:eq+n_eq-1, vbase2+1:vbase2+problem.state_dof) = j2;
        end
        if nargout > 2
            W(eq:eq+n_eq-1,eq:eq+n_eq-1) = w;
        end
        eq = eq+n_eq;
    end


    F = F';

