function [f,w,jc,j] = meas_depth_norot(common,state,obs)
    % Careful, X in NED
    f = state.X(3) + obs.depth;
    w = 1.0/0.25^2;
    if nargout > 2
        jc = zeros(1,3);
        j = zeros(1,4); 
        j(3) = 1;
    end


