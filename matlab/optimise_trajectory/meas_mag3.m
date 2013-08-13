function [f,w,jc,j] = meas_mag3(common,state,obs)
    % Ideal, we would write it this way
    % f = rotz(pi/2)*state.R' * obs.Bref - obs.B;
    % However, it is nicer to rewrite it as below (to express the jacobian) and
    % it does not change the behaviour of the error
    f = state.R*obs.Bref - rotz(pi/2)*obs.B ;
    if nargout > 1
        w = eye(3) * 10.0;
    end

    if nargout > 2
        j = zeros(3,7);
        jc = [];
        j(1:3,5:7) = -cross_mat(state.R*obs.Bref);
    end


