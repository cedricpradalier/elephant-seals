clear all
close all
addpath('igrf');
addpath('rotations');
addpath('optimisation');

load('data/gps.mat');
load('data/preload.mat');

N = size(preload,1);
ts = preload(1:N,1) + preload(1:N,2);
A = preload(1:N,3:5);
M = preload(1:N,6:8);
depth = preload(1:N,9);
RPY = preload(1:N,15:17);

t = (ts - ts(1))*24*3600;


common = {}
common.Bscale = 100.0*ones(3,1);
common.k_depth = 2.0 / 600.;
Bned = igrf(ts(1),gps(1,3),gps(1,2),0);
% Benu = [Bned(2), Bned(1), -Bned(3)];
common.Bref = Bned';
common.Aref = [0;0;-10];
state = {}
obs = {}
for i=1:N
    state{i}.ts = ts(i);
    state{i}.depth = depth(i);
    state{i}.P = 0.0;
    state{i}.R = rpy(RPY(i,1),RPY(i,2),0);
    obs{i}.A = A(i,:)';
    obs{i}.B = M(i,:)';
end

problem=struct();
problem.common_dof = 4;
problem.state_dof = 4;
problem.f_state={@meas_accel_cal,@meas_mag_cal};
% problem.f_state={};
problem.f_delta={@cont_propulsion_cal};
problem.f_delta={};
problem.f_async={};
problem.f_async2={};

mu = 1;
scale = 1.1;
for iter=1:1000
    disp 'Building optimisation problem'
    [E,W,J]=cost(common,state,obs,problem);
    current_error = E'*W*E

    while mu < 20
        disp 'Solving for delta'
        delta = -(J'*W*J + mu*eye(size(J,2))) \ (J'*W*E);
        if (norm(delta/size(delta,1))<1e-3) || (mu >= 20)
            break
        end

        disp 'Computing new state'
        dcommon = {};
        dstate = {};
        new_state = state;
        new_common = common;
        dcommon.Bscale = delta(1:3);
        dcommon.k_depth = delta(4);
        new_common.Bscale = common.Bscale + dcommon.Bscale;
        new_common.k_depth = common.k_depth + dcommon.k_depth;
        for i=1:N
            vbase = problem.common_dof + (i-1)*problem.state_dof;
            dstate{i}.R = exp_mat(delta(vbase+1:vbase+3));
            dstate{i}.P = delta(vbase+4);

            new_state{i}.R = dstate{i}.R * state{i}.R ;
            new_state{i}.P = state{i}.P + dstate{i}.P ;
        end
        
        Ap = zeros(N,3); Bp = zeros(N,3);
        for i=1:N
            Ap(i,:) = (new_state{i}.R * obs{i}.A - [0;0;new_common.k_depth]*new_state{i}.depth)';
            Bp(i,:) = (new_state{i}.R * diag(common.Bscale) * obs{i}.B )';
        end
        figure(1); plot(t,Ap(:,1),'b',t,Ap(:,2),'r',t,Ap(:,3),'g');
        figure(2); plot(t,Bp(:,1),'b',t,Bp(:,2),'r',t,Bp(:,3),'g');
        drawnow();
        % pause

        [E,W,J]=cost(new_common,new_state,obs,problem);
        new_error = E'*W*E
        if (mu==0) || (new_error < current_error)
            state = new_state;
            common = new_common;
            mu = mu / scale
            disp 'accepted'
            break
        else
            disp 'increasing mu'
            mu = mu * scale
        end
    end

    % state = new_state;
    % common = new_common;
    norm(delta/size(delta,1))
    if (norm(delta/size(delta,1))<1e-3) || (mu >= 20)
        break
    end
end

RPY=zeros(N,3);
P=zeros(N,1);
for i=1:N
    P(i) = state{i}.P;
    RPY(i,:) =  SpinCalc('DCMtoEA321',state{i}.R,1,0) ;
end

% save -mat data/state.mat state common obs
