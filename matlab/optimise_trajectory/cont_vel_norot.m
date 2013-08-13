function [f,w,jc,j] = cont_vel_norot(common,state1,state2,obs1,obs2)
    f = state2.V - state1.V;
    if nargout > 1
        w = 1/0.2^2;
    end
    if nargout > 2
        jc = zeros(1,3);
        j = zeros(1,8);
        j(4) = -state1.V; j(8) = state2.V;
    end


