function [f,w,jc,j,s] = meas_gps_norot(gps,common,states)
    s = gps(4);
    if s > size(states,2)
        f = []
    else 
        X = states{gps(4)}.X(1:2);
        f = [X(1) - gps(7); X(2) - gps(6)];
    end
    if nargout > 1
        w = 1e6 * eye(2);
    end
    if nargout > 2
        jc = zeros(2,3);
        j = zeros(2,4);
        j(1,1) = 1;
        j(2,2) = 1;
    end


