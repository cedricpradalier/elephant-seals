function [f,w,jc,j] = cont_propulsion_cal(common,state1,state2,obs1,obs2)
    f = state2.P - state1.P;
    if nargout > 1
        w = 1e3 * eye(1);
    end
    if nargout > 2
        jc = zeros(1,4);
        j = zeros(1,8);
        j(1, 8) = 1;
        j(1, 4) = -1;
    end


