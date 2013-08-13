function R = exp_mat(x)
    cross = cross_mat(x);
    theta = norm(x);
    if (abs(theta)<1e-4)
        R = eye(3);
    else
        R = eye(3) + cross*sin(theta)/theta + cross*cross*(1-cos(theta))/theta^2;
    end
