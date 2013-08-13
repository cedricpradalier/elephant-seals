function [f,w,jc,j] = meas_mag_cal(common,state,obs)
    f = common.Bscale .* obs.B - state.R'*common.Bref ;
    if nargout > 1
        w = eye(3) * 0.10;
    end

    if nargout > 2
        jc = zeros(3,4);
        jc(1:3,1:3) = diag(obs.B);
        j = zeros(3,4);
        j(1:3,1:3) = -state.R'*cross_mat(common.Bref);
    end


