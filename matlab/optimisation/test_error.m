
theta=pi/3;
R1 = [cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0; 0 0 1];
theta=pi/4;
R2 = [cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta) ];
R = R1*R2;
beta = (trace(R) - 1)/2;
f_0 = acos(beta)
acos_prime=@(x)-1/sqrt(1-x^2);

sigma = -[R(3,2)-R(2,3), R(1,3)-R(3,1), R(2,1)-R(1,2)]/2;
J = acos_prime(beta) * sigma;

B = J'*f_0;
A = J'*J;
rank(A)

% delta = -inv(J'*J)*J*f_0
% delta = -inv(A)*B

[U, S, V] = svd(A);

C = -U'*B;
n = rank(S);
delta_u = inv(S(1:n,1:n)) * C(1:n);
delta_v = [delta_u;zeros(3-n,1)];
delta = V*delta_v


Rup = exp_mat(delta)*R



