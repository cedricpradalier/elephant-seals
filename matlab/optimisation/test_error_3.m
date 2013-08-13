
theta=pi/3;
R = [cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0; 0 0 1];

Ri = R
for i=1:2
    V = log_mat(Ri);
    if norm(V) < 1e-4
        break
    end

    F_0 = V;
    J = Slog_inv(V);
    delta = -inv(J'*J)*J'*F_0;
    Ri = exp_mat(delta)*Ri

end


