function X = lm_basic(f,X0, epsilon)
    Xprev = X0 + ones(size(X0));
    X = X0;
    mu = 1;
    scale = 1.1;
    while norm(X - Xprev) > epsilon
        [F,J] = f(X);
        % F

        % solve for delta
        delta = -inv(J'*J + mu*eye(length(J))) * J' * F;
        Xcandidate = X + delta;
        value = f(Xcandidate);
        if (value < F)
            Xprev = X;
            X = Xcandidate;
            mu = mu / scale;
        else 
            mu = mu * scale;
        end
    end
