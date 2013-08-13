function Sinv = Slog_inv(v)
    phi = norm(v);
    v_x = cross_mat(v);
    if phi < 1e-4
        Sinv = eye(3) + v_x/2 ;
    else
        a = (1 - phi*cot(phi/2)/2)/phi^2;
        Sinv = eye(3) + v_x/2 + a * v_x * v_x;
    end
