lambda = 0.01;
P = [0 2; 2 0; 4 -1 ; 6 2];
Preal = P + randn(4,2)*0.1;

l = [norm(P(1,:)-P(2,:))  norm(P(2,:)-P(3,:))  norm(P(3,:)-P(4,:))];
lreal = l + randn(1,3)*0.1;

Z = [lreal   Preal(1,:) Preal(4,:) P(:,2)'];

lW = 2.00;
pW = 1.00;
dW = 0.01;
W = [ lW*ones(1,3) pW*ones(1,4) dW*ones(1,4)];

X = [zeros(4,1) randn(4,1)];
X = reshape(X',1,8);
ds = ones(1,8);
while norm(ds) > 1e-5

    Z
    J = jac(X);
    O = obs(X);
    E = (Z - O) ./ W;
    dX = (J\(E'))';
    ds = norm(dX)

    X = X + (lambda * dX)
    % pause
end

