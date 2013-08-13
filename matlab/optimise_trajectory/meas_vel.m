function [f,w,jc,j] = meas_vel(common,state,obs)
    if obs.dive == 0
        f = state.V;
        w = 1/0.5^2;
    elseif obs.has_vel
        f = state.V - obs.vel;
        w = 2.0;
    else 
        f = [];
        w = 1.0;
    end

    if nargout > 2
        jc = [];
        j = zeros(1,7);
        j(4) = 1;
    end

