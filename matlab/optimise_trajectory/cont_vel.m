function [f,w,jc,j] = cont_vel(common,state1,state2,obs1,obs2)
    f = state2.V - state1.V;
    if nargout > 1
        w = 1/0.2^2;
    end
    if nargout > 2
        jc = [];
        j = zeros(1,14);
        j(4) = -1; j(11) = 1;
    end


