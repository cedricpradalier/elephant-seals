function [H,J] = potential_angle(X,Z)

a=1.0;
b=1.0;
c=1.0;
s=sum([a,b,c]');
a=a/s;
b=b/s;
c=c/s;

n = length(X)/7;
x = {}
for i=0:n-1
    istart=4*i+1
    q_est = X(istart:istart+4)
    q_meas = Z(istart:istart+4)
    istart=4*n+1+3*i
    w_est = X(istart:istart+3)
    x{i+1}=struct("q_est",q_est,"w_est",w_est,"q_meas",q_meas)
    x{i+1}.err_z = quat_mat(quat_mul(q_est,quat_inv(q_meas)))
end
for i=1:n-1
    x{i}.err_m = quat_mat(quat_mul(quat_mul(x{i}.w_est,x{i}.q_est),   ...
        quat_inv(x{i+1}.q_est)))
    x{i}.err_w = x{i}.w_est - x{i+1}.w_est
end

Iz = 0
Im = 0
Iw = 0

for i=1:n
    Iz = Iz + acos((tr(x{i}.err_z)-1)/2)^2
end
for i=1:n-1
    Im = Im + acos((tr(x{i}.err_m)-1)/2)^2
    Iw = Iw + norm(x{i}.err_w)^2
end

H = a*Iz + b*Im + c*Iw;
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


