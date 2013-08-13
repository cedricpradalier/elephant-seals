% clear all
close all
addpath('../..');

Bref = [ -0.29446  0.19735  0.93506]';

xmax = 100.0;
ymax = 50.0;
zmax = 200.0;

N = 300;

ts = linspace(0,N,N)';
x = ts * xmax / N;
y = sin(ts * 2 * pi / N).^3 * ymax;
z = 0.5 * (cos(ts * 2 * pi / N) - 1) * zmax;

X = [x y z];
dX = ([diff(X) ; 0 0 0] + [0 0 0; diff(X)])/2;
V = sqrt(sum((dX .* dX)'))';

has_vel = ones(N,1);
vel = V + randn(N,1)*0.2;
depth = z + randn(N,1)*0.1;

roll = zeros(N,1);
pitch = atan2(-dX(:,3),hypot(dX(:,1),dX(:,2)));
yaw = atan2(dX(:,2),dX(:,1));

An = zeros(N,3);
Mn = zeros(N,3);
vdt = zeros(N,3);
rolln = roll + randn(N,1)*(3*pi/180);
pitchn = pitch + randn(N,1)*(3*pi/180);
yawn = yaw + randn(N,1)*(3*pi/180);

for i=1:N
    R = rpy(rolln(i),pitchn(i),yawn(i)+0.2);
    vdt(i,:) = (R * [vel(i);0;0])'; % dt = 1
    R = rpy(roll(i),pitch(i),yaw(i));
    An(i,:) = (R' * [0;0;-1])' + randn(1,3)*0.1;
    An(i,:) = An(i,:) / norm(An(i,:));
    Mn(i,:) = (R' * Bref)' + randn(1,3)*0.1;
    Mn(i,:) = Mn(i,:) / norm(Mn(i,:));
end

Xint = [0 0 0;cumsum(vdt(1:end-1,:))];

dive = abs(depth) > 1.0;

ts = ts / (24*3600);


gps = [0 71.453 -49.184 1 0.0 0.0 0.0 43 Bref(2) Bref(1) -Bref(3); 
ts(end) 71.454 -49.183 N 0.0 X(end,1) X(end,2) 43 Bref(2) Bref(1) -Bref(3)];
save -mat gps_sim.mat gps


preload = [ts, An, Mn, depth, vel, has_vel, Xint, rolln, pitchn, yawn, vel, dive, X, V, roll, pitch, yaw];
preload_titles ={};
preload_titles{1} = 'datenum';
preload_titles{2} = 'Ax';
preload_titles{3} = 'Ay';
preload_titles{4} = 'Az';
preload_titles{5} = 'Mx';
preload_titles{6} = 'My';
preload_titles{7} = 'Mz';
preload_titles{8} = 'depth';
preload_titles{9} = 'vel';
preload_titles{10} = 'has_vel';
preload_titles{11} = 'X';
preload_titles{12} = 'Y';
preload_titles{13} = 'Z';
preload_titles{14} = 'roll';
preload_titles{15} = 'pitch';
preload_titles{16} = 'yaw';
preload_titles{17} = 'V';   
preload_titles{18} = 'dive_status';   
save -mat preload_sim.mat preload preload_titles


