function [f,w,jc,j] = kinematic(common,state1,state2,obs1,obs2)
    dt = (state2.ts-state1.ts) * (24*3600);
    V = [state1.V*dt;0;0];
    Rp = state1.R;
    f = state2.X - (state1.X + Rp*V);
    if nargout > 1
        % w = inv(Rp' * diag([1.0^2, 0.5^2, 0.5^2]) * Rp);
        w = 1e4*eye(3);
    end
    if nargout > 2
        jc = [];
        j = zeros(3,14);
        j(1:3, 8:10) = eye(3);
        j(1:3, 1:3) = -eye(3);
        j(1:3, 4) = -Rp(:,1)*dt;
        j(1:3, 5:7) = cross_mat(Rp * V);
    end

