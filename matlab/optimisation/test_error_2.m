
theta1=pi/3;
R1 = [cos(theta1) -sin(theta1) 0;sin(theta1) cos(theta1) 0; 0 0 1];
theta2=pi/6;
R2 = [cos(theta2) 0 sin(theta2);0 1 0; -sin(theta2) 0 cos(theta2)];
dR=R2*R1';

R1_obs = exp_mat(randn(3,1)*1e-2)*R1;
R2_obs = exp_mat(randn(3,1)*1e-2)*R2;
dR_obs = exp_mat(randn(3,1)*1e-2)*dR;

R1_i = R1_obs;
R2_i = R2_obs;

mu = 1;
scale = 1.1;
for i=1:1000

    E1 = R1_i*R1_obs';
    E2 = R2_i*R2_obs';
    Ed = (R2_i * R1_i') * dR_obs';
    Edback = (R1_i * R2_i') * dR_obs';

    vE1 = log_mat(E1);
    vE2 = log_mat(E2);
    vEd = log_mat(Ed);

    F_0 = [vE1;vE2;vEd];
    cost = sum(F_0.*F_0)
    J = zeros(9,6);
    J(1:3,1:3) = Slog_inv(E1);
    J(4:6,4:6) = Slog_inv(E2);
    J(7:9,1:3) = -dR_obs*Slog_inv(Edback);
    J(7:9,4:6) = Slog_inv(Ed);

    delta = -inv(J'*J + mu*eye(size(J,2)))*J'*F_0;
    if norm(delta)<1e-6
        break
    end
    cR1_i = exp_mat(delta(1:3))*R1_i;
    cR2_i = exp_mat(delta(4:6))*R2_i;
    RR1_i = [R1 cR1_i]
    RR2_i = [R2 cR2_i]
    
    E1 = cR1_i*R1_obs';
    E2 = cR2_i*R2_obs';
    Ed = (cR2_i * cR1_i') * dR_obs';
    Edback = (cR1_i * cR2_i') * dR_obs';

    vE1 = log_mat(E1);
    vE2 = log_mat(E2);
    vEd = log_mat(Ed);

    F_0 = [vE1;vE2;vEd];
    new_cost = sum(F_0.*F_0);
    if new_cost < cost
        R1_i = cR1_i;
        R2_i = cR2_i;
        mu = mu / scale
    else
        mu = mu * scale
    end

end




