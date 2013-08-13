function [f,w,jc,j] = meas_accel3(common,state,obs)
    f = state.R * obs.A - [0;0;-1];
    if nargout > 1
        w = eye(3) * 10.0;
    end

    if nargout > 2
        j = zeros(3,7); jc = [];
        j(1:3,5:7) = -cross_mat(state.R*obs.A);
    end


