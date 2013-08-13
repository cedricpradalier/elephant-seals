function X = gauss_newton(f,X0, epsilon)
    Xprev = X0 + ones(size(X0))
    X = X0
    while norm(X - Xprev) > epsilon
        [F,J] = f(X);
        F
        % solve for delta
        delta = -inv(J'*J) * J' * F;
        Xprev = X;
        X = X + 0.01*delta';
    end
