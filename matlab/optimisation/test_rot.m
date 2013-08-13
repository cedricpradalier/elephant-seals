t = [0:36]';
w = pi/18
R = 1.0
n = length(t)


x = R * cos(w * t);
y = R * sin(w * t);
theta = w * t + pi/2;
theta_ref = theta;
w_ref = w * ones(length(t),1);

figure(1)
plot(x,y,'-b')
hold on
quiver(x,y,cos(theta),sin(theta),0)
hold off
axis equal


theta_measure = theta + randn(length(t),1)*pi/(4*18);
figure(2)
plot(x,y,'-b')
hold on
quiver(x,y,cos(theta_measure),sin(theta_measure),0)
hold off
axis equal


theta_est = theta_measure;
w_est = [diff(theta_measure);0];

figure(3);
plot(t,w_ref,'b',t,w_est,'r')

figure(4);
plot(t,theta_est - theta_ref,'b')

X = [theta_est;zeros(size(t))];
Z = theta_measure;

f=@(x)potential(x,Z);
% matlab only
% Xopt = lsqnonlin(f,X,[],[],OPTIONS);

Xopt = lm_basic(f,X,1e-4)


theta_opt=Xopt(1:n)
w_opt=Xopt(n+1:2*n)

figure(3);
plot(t,w_ref,'b',t,w_est,'r',t,w_opt,'g')

figure(4);
plot(t,theta_est - theta_ref,'b',t,theta_opt-theta_ref,'g')

figure(1)
plot(x,y,'-b')
hold on
quiver(x,y,cos(theta_opt),sin(theta_opt),0)
hold off
axis equal

