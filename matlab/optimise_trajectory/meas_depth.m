function [f,w,jc,j] = meas_depth(common,state,obs)
    if obs.dive == 0
        f = state.X(3) - obs.depth;
        w = 1.0;
    else
        f = state.X(3) - obs.depth;
        w = 1.0/0.25^2;
    end
    if nargout > 2
        j = zeros(1,7); jc=[];
        j(3) = 1;
    end


