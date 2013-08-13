function [f,w,jc,j] = meas_vel_norot(common,state,obs)
    if abs(state.X(3)) < 2.0
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
        jc = zeros(1,3);
        j = zeros(1,4);
        j(4) = state.V;
    end

