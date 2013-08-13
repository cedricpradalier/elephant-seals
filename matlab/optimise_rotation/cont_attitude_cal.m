function [f,w,jc,j] = cont_attitude_cal(common,state1,state2,obs1,obs2)
    Ed = (state2.R * state1.R');
    f = (log_mat(Ed));
    if nargout > 1
        w = 1e1 * eye(3);
    end
    if nargout > 2
        jc = zeros(3,4);
        j = zeros(3,8);
        j(1:3, 5:7) = Slog_inv(Ed);
        j(1:3, 1:3) = -Slog_inv(Ed');
    end


