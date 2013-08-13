function [H,J] = potential(X,Z)

a=1.0;
b=1.0;
c=1.0;
s=sum([a,b,c]');
a=a/s;
b=b/s;
c=c/s;

n = length(X)/2;
theta_est = X(1:n);
w_est = X(n+1:2*n);
theta_meas = Z;

Iz = theta_est - theta_meas;
Im = theta_est(1:n-1) + w_est(1:n-1) - theta_est(2:n);
Iw = w_est(1:n-1) - w_est(2:n);

H = a*sum(Iz .* Iz) + b*sum(Im .* Im) + c*sum(Iw .* Iw);
H = H / 2.0;

if nargout > 1   % Two output arguments
    Jt=zeros(n,1);
    Jw=zeros(n,1);
    % First the jacobian of H wrt theta
    Jt(1) = a*(theta_est(1) - theta_meas(1)) ...
            + b*(theta_est(1) + w_est(1) - theta_est(2));
    for i=2:n-1
        Jt(i) = a*(theta_est(i) - theta_meas(i)) ...
            + b*(theta_est(i) + w_est(i) - theta_est(i+1))...
            - b*(theta_est(i-1) + w_est(i-1) - theta_est(i));
    end
    Jt(n) = a*(theta_est(n) - theta_meas(n)) ...
            - b*(theta_est(n-1) + w_est(n-1) - theta_est(n));

    % Then the jacobian wrt w
    Jw(1) = b*(theta_est(1) + w_est(1) - theta_est(2)) ...
            - c*(w_est(2) - w_est(1));
    for i=2:n-1
        Jw(i) = b*(theta_est(i) + w_est(i) - theta_est(i+1)) ...
            + c*(w_est(i)-w_est(i-1)) - c*(w_est(i+1)-w_est(i));
    end
    Jw(n) = c*(w_est(n) - w_est(n-1));

    J = [Jt;Jw]';
end


