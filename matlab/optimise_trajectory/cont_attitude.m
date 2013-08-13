function [f,w,jc,j] = cont_attitude(common,state1,state2,obs1,obs2)
    Ed = (state2.R * state1.R');
    f = (log_mat(Ed));
    if nargout > 1
        w = 1e1 * eye(3);
    end
    if nargout > 2
        jc = [];
        j = zeros(3,14);
        j(1:3, 12:14) = Slog_inv(Ed);
        j(1:3, 5:7) = -Slog_inv(Ed');
    end


