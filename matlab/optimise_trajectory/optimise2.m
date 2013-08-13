% clear all
% close all
addpath('../..');
addpath('../../optimisation');
acos_prime=@(x)-1/sqrt(1-x^2);

load('gps.mat');

load('preload.mat');

N = 500; % size(preload,1);
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


X(:,1:2) = X(:,1:2) + ones(size(X,1),1)*gps(1,6:7);
common = {}
state = {}
obs = {}
for i=1:N
    state{i}.ts = ts(i);
    state{i}.X = X(i,:)';
    state{i}.V = V(i,:);
    state{i}.R = rpy(RPY(i,1),RPY(i,2),RPY(i,3));
    obs{i}.A = An(i,:)';
    obs{i}.B = Mn(i,:)';
    Bned = gps(1,9:11);
    Benu = [Bned(2), Bned(1), -Bned(3)];
    obs{i}.Bref = Bned';
    obs{i}.dive = dive(i);
    obs{i}.vel = vel(i);
    obs{i}.depth = depth(i);
    obs{i}.has_vel = has_vel(i);
end

problem=struct();
problem.state_dof = 7;
problem.common_dof = 0;
problem.f_state={@meas_depth,@meas_accel3,@meas_mag3,@meas_vel};
% problem.f_state={};
problem.f_delta={@cont_vel,@cont_attitude,@kinematic};
% problem.f_delta={@kinematic};
% problem.f_delta={};
problem.f_async={};
for i=1:Ngps
    problem.f_async{end+1}=@(x)meas_gps(gps(i,:),x);
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
        new_state = state;
        for i=1:N
            vbase = (i-1)*problem.state_dof;
            dstate{i}.X = delta(vbase+1:vbase+3);
            dstate{i}.V = delta(vbase+4);
            dstate{i}.R = exp_mat(delta(vbase+5:vbase+7));

            new_state{i}.X = state{i}.X + dstate{i}.X ;
            new_state{i}.V = state{i}.V + dstate{i}.V ;
            new_state{i}.R = dstate{i}.R * state{i}.R ;
        end
        X1 = zeros(N,8); B1 = zeros(N,3);
        X2 = zeros(N,8); B2 = zeros(N,3);
        for i=1:N
            X1(i,1:4) = [state{i}.ts,state{i}.X'];
            X1(i,5) = state{i}.V;
            % X1(i,6:8) = mat2rpy(state{i}.R)
            B1(i,:) = (state{i}.R * obs{i}.A)';
            X2(i,1:4) = [state{i}.ts,new_state{i}.X'];
            X2(i,5) = new_state{i}.V;
            % X2(i,6:8) = mat2rpy(new_state{i}.R)
            B2(i,:) = (new_state{i}.R * obs{i}.A)';
        end
        plot(X1(:,2),X1(:,3),'b',X2(:,2),X2(:,3),'r');
        drawnow();

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

    if (norm(delta/size(delta,1))<1e-3) || (mu >= 20)
        break
    end
end

X1 = zeros(N,8); B1 = zeros(N,3);
X2 = zeros(N,8); B2 = zeros(N,3);
for i=1:N
    X1(i,1:4) = [state{i}.ts,state{i}.X'];
    X1(i,5) = state{i}.V;
    % X1(i,6:8) = mat2rpy(state{i}.R)
    B1(i,:) = (state{i}.R * obs{i}.A)';
    X2(i,1:4) = [state{i}.ts,new_state{i}.X'];
    X2(i,5) = new_state{i}.V;
    % X2(i,6:8) = mat2rpy(new_state{i}.R)
    B2(i,:) = (new_state{i}.R * obs{i}.A)';
end

plot(X(:,1),X(:,2),'b',X2(:,2),X2(:,3),'r');
hold on
S=1:10:N;
quiver(X1(S,2),X1(S,3),B1(S,1),B1(S,2),'b');
quiver(X2(S,2),X2(S,3),B2(S,1),B2(S,2),'r');
hold off

