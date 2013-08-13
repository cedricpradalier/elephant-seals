function [f,w,jc,j] = meas_accel_cal(common,state,obs)
    f = state.R * (obs.A - [state.P;0;0]) - (common.Aref + [0;0;common.k_depth * state.depth]);
    if nargout > 1
        w = eye(3) * 10.0;
    end

    if nargout > 2
        jc = zeros(3,4);
        jc(3,4) = -state.depth;
        j = zeros(3,4); 
        j(1:3,1:3) = -cross_mat(state.R*obs.A);
        j(1:3,4) = -state.R(:,1);
    end


