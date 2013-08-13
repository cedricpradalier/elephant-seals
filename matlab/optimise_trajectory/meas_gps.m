function [f,w,jc,j,s] = meas_gps(gps,common,states)
    s = gps(4);
    if s > size(states,2)
        f = []
    else 
        f = states{gps(4)}.X(1:2) - gps(6:7)';
    end
    if nargout > 1
        w = 1e6 * eye(2);
    end
    if nargout > 2
        jc = [];
        j = zeros(2,7);
        j(1,1) = 1;
        j(2,2) = 1;
    end


