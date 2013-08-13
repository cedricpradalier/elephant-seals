% clear all
% close all
addpath('../..');
addpath('../../optimisation');

% This optimisation assumes that optimise3 has been run to compute the proper
% rotation matrices



load('gps.mat');

load('state.mat');

load('preload.mat');

N = size(preload,1);
ts = preload(1:N,1);
An = preload(1:N,2:4);
depth = preload(1:N,8);
vel = preload(1:N,9);
has_vel = preload(1:N,10);
V = preload(1:N,17);
dive = preload(1:N,18);

% Caution, R is tuned in NED frame, so we stay there
X = zeros(N,3);
for i=2:N
    dt = (ts(i)-ts(i-1))*24*3600;
    X(i,:) = (X(i-1,:)'+ state{i}.R * [V(i-1) * dt;0;0])';
end
X(:,3) = -(depth - max(depth));
dive(find(abs(X(:,3))<2.0)) = 0;
V(find(dive==0)) = 0.0;
% X=[X(:,2),X(:,1),depth];

disp 'Preloaded and computed trajectory'

% Commented to work in simulation
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
Ngps = size(gps,1);


X(:,1) = X(:,1) + gps(1,7);
X(:,2) = X(:,2) + gps(1,6);
disp 'ready to build optimisation problem'

% Rotation between IMU attitude and animal displacement
common.R_attach = eye(3);
for i=1:N
    state{i}.ts = ts(i);
    state{i}.X = X(i,:)';
    state{i}.V = V(i,:);
    obs{i}.dive = dive(i);
    obs{i}.vel = vel(i);
    obs{i}.depth = depth(i);
    obs{i}.has_vel = has_vel(i);
end

problem=struct();
problem.state_dof = 4;
problem.common_dof = 3;
problem.f_state={@meas_depth_norot,@meas_vel_norot};
% problem.f_state={};
problem.f_delta={@cont_vel_norot,@kinematic_norot};
% problem.f_delta={@kinematic};
% problem.f_delta={};
problem.f_async={};
for i=1:Ngps
    problem.f_async{end+1}=@(c,x)meas_gps_norot(gps(i,:),c,x);
end
problem.f_async2={};

mu = 0;
scale = 1.1;
for iter=1:1000
    disp 'Building optimisation problem'
    [E,W,J]=cost(common,state,obs,problem);
    current_error = E'*W*E

    % while mu < 20
        disp 'Solving for delta'
        delta = -(J'*W*J + mu*eye(size(J,2))) \ (J'*W*E);
        if (norm(delta/size(delta,1))<1e-6) || (mu >= 20)
            break
        end

        disp 'Computing new state'
        dstate = {};
        dcommon = {};
        new_state = state;
        new_common = common;
        for i=1:N
            vbase = problem.common_dof + (i-1)*problem.state_dof;
            dstate{i}.X = delta(vbase+1:vbase+3);
            dstate{i}.V = delta(vbase+4);
            dcommon.R = exp_mat(delta(1:3));

            new_state{i}.X = state{i}.X + dstate{i}.X ;
            new_state{i}.V = exp(log(max(state{i}.V,1e-5)) + dstate{i}.V) ;
            new_common.R_attach = dcommon.R * common.R_attach;
        end


        % [E,W,J]=cost(new_state,obs,problem);
        % new_error = E'*W*E
        % if (mu==0) || (new_error < current_error)
        %     state = new_state;
        %     mu = mu / scale
        %     disp 'accepted'
        %     break
        % else
        %     disp 'increasing mu'
        %     mu = mu * scale
        % end
    % end
    state = new_state;
    common = new_common;
    Y = zeros(N,4);
    for i=1:N
        Y(i,1:3) = state{i}.X';
        Y(i,4) = state{i}.V;
    end
    plot(Y(:,2),Y(:,1),'b',X(:,2),X(:,1),'r',gps(:,6),gps(:,7),'*g','markersize',10);
    drawnow();

    if (norm(delta/size(delta,1))<1e-3) || (mu >= 20)
        break
    end
end

