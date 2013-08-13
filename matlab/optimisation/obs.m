function O = obs(X)
    O = [norm(X(1:2)-X(3:4))  norm(X(3:4)-X(5:6))  norm(X(5:6)-X(7:8))  X(1:2)  X(7:8) X(2:2:8)];

