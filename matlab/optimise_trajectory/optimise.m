clear all
close all
addpath('../..');
addpath('../../optimisation');
acos_prime=@(x)-1/sqrt(1-x^2);

load('gps.mat');

load('preload.mat');

N = 6; % size(preload,1);
ts = preload(:,1);
An = preload(:,2:4);
Mn = preload(:,5:7);
depth = preload(:,8);
vel = preload(:,9);
has_vel = preload(:,10);
X = preload(:,11:13);
RPY = preload(:,14:16);
V = preload(:,17);
dive = preload(:,18);



ts_surf=ts(find(dive==0),:);
for i=1:size(gps,1)
    [v,w] = min(abs(ts_surf-gps(i,1)));
    w = find(ts==ts_surf(w));
    gps(i,5) = v;
    if abs(v)<(30./(24*60))
        gps(i,4) = w;
    else
        gps(i,4) = 0;
    end
end
gps = gps(find(gps(:,4)>0),:);
Ngps = 1 % size(gps,1);


X(:,1:2) = X(:,1:2) + ones(size(X,1),1)*gps(1,6:7);
state = {}
obs = {}
for i=1:N
    state{i}.X = X(i,:);
    state{i}.V = V(i,:);
    state{i}.R = rpy(RPY(i,1),RPY(i,2),0);
    obs{i}.A = An(i,:);
    obs{i}.B = -Mn(i,:);
    Bned = gps(1,9:11);
    obs{i}.Bref = [Bned(2), Bned(1), -Bned(3)];
    obs{i}.dive = dive(i);
    obs{i}.vel = vel(i);
    obs{i}.has_vel = has_vel(i)
end

N_var = 7
N_eq = 11
Ngps_eq = 2

E = zeros(N_eq*N + Ngps_eq*Ngps,1);
J = sparse(N_eq*N + Ngps_eq*Ngps,N_var*N);
for i=1:N
    cbase=N_eq*(i-1);
    vbase=N_var*(i-1);
    vbase_prev = vbase - N_var;
    % Depth
    if dive(i) == 0
        E(cbase+1) = state{i}.X(3);
    else
        E(cbase+1) = state{i}.X(3) - depth(i);
    end
    J(cbase+1, vbase+3) = 1;
    % Accelerometer
    dot = An(i,:) * state{i}.R' * [0;0;-1];
    E(cbase+2) = acos(dot);
    if (abs(dot-1)>1e-4)
        jac = - acos_prime(dot) * cross_mat(An(i,:)) * state{i}.R'*[0;0;-1];
    else
        jac = zeros(3,1);
    end
    J(cbase+2, vbase+5) = jac(1);
    J(cbase+2, vbase+6) = jac(2);
    J(cbase+2, vbase+7) = jac(3);
    
    % Magnetometer
    Bned = gps(1,9:11);
    Benu = [Bned(2), Bned(1), -Bned(3)];
    dot = -Mn(i,:) * state{i}.R' * Benu';
    if (abs(dot-1)>1e-4)
        jac = - acos_prime(dot) * cross_mat(Mn(i,:)) * state{i}.R'*Benu';
    else
        jac = zeros(3,1);
    end
    E(cbase+3) = acos(dot);
    J(cbase+3, vbase+5) = jac(1);
    J(cbase+3, vbase+6) = jac(2);
    J(cbase+3, vbase+7) = jac(3);
    
    % Velocity if available
    if dive(i) == 0
        E(cbase+4) = state{i}.V;
        J(cbase+4, vbase+4) = 1;
    elseif has_vel(i)
        E(cbase+4) = state{i}.V - vel(i);
        J(cbase+4, vbase+4) = 1;
    end

    if i>1
        % Rotation continuity
        E(cbase+5:cbase+7) = (log_mat(state{i}.R * state{i-1}.R'));
        Ed = (state{i}.R * state{i-1}.R');
        J(cbase+5:cbase+7, vbase+5:vbase+7) = Slog_inv(Ed);
        J(cbase+5:cbase+7, vbase_prev+5:vbase_prev+7) = -Slog_inv(Ed');

        % Velocity continuity
        E(cbase+8) = state{i}.V - state{i-1}.V;
        J(cbase+8, vbase+4) = 1;
        J(cbase+8, vbase_prev+4) = -1;

        % % Kinematic model
        dt = ts(i)-ts(i-1);
        V = [state{i-1}.V*dt;0;0];
        Rp = state{i-1}.R;
        E(cbase+9:cbase+11) = state{i}.X' - (state{i-1}.X' + Rp*V);
        J(cbase+9:cbase+11, vbase+1:vbase+3) = eye(3);
        J(cbase+9:cbase+11, vbase_prev+1:vbase_prev+3) = -eye(3);
        J(cbase+9:cbase+11, vbase_prev+4) = -Rp(:,1);
        J(cbase+9:cbase+11, vbase_prev+5:vbase_prev+7) = -Rp*cross_mat(V)*Rp';
    end
end

for i=1:Ngps
    % GPS anchor points
    cbase = N_eq*N + (i-1)*Ngps_eq;
    E(cbase+1:cbase+2) = state{gps(i,4)}.X(1:2)' - gps(i,6:7)';
    J(cbase+1,N_var*(gps(i,4)-1)+1) = 1;
    J(cbase+2,N_var*(gps(i,4)-1)+2) = 1;
end


delta = J \ E;

dstate = {};
new_state = {}
for i=1:N
    vbase = (i-1)*N_var;
    dstate{i}.X = delta(vbase+1:vbase+3)';
    dstate{i}.V = delta(vbase+4);
    dstate{i}.R = exp_mat(delta(vbase+5:vbase+7));

    new_state{i}.X = state{i}.X + dstate{i}.X ;
    new_state{i}.V = state{i}.V + dstate{i}.V ;
    new_state{i}.R = dstate{i}.R * state{i}.R ;
end

