clear all
close all
addpath('igrf');

load testorientation.mat
data = acou;
B = igrf(time(1),46+7/60,0+20/60,0);
B = B / (sqrt(sum(B.*B)))

N = size(data,1)
M = data(:,4:6);
A = data(:,7:9);
nA=sqrt(sum((A.*A)')');
nM=sqrt(sum((M.*M)')');
An = A .* ((1./nA) * [1 1 1]);
Mn = M .* ((1./nM) * [1 1 1]);

pitch=asin(An(:,1));
Aroll = zeros(N,3);
for i=1:N
    R = rpy(0,-pitch(i),0);
    Aroll(i,:) = (R * An(i,:)')';
end
roll=-asin(Aroll(:,2));
Anone = zeros(N,3);
Anone2 = zeros(N,3);
Mr = zeros(N,3);
for i=1:N
    R = rpy(roll(i),0,0);
    Anone(i,:) = (R * Aroll(i,:)')';
    R = rpy(roll(i),pitch(i),0);
    Anone2(i,:) = (R * An(i,:)')';
    Mr(i,:) = (R * Mn(i,:)')';
end
yaw=fmod(atan2(Mr(:,2),Mr(:,1)),2*pi);
Mnone = zeros(N,3);
for i=1:N
    % R = rpy(0,0,-yaw(i));
    R = rpy(roll(i),pitch(i),yaw(i));
    Mnone(i,:) = (R * Mr(i,:)')';
end

B = igrf(time(1),46+7/60,0+20/60,0);
B = B / norm(B);
Bmat = inv([B(1,1) -B(1,2);B(1,2) B(1,1)]);
csy=(Bmat * Mr(:,1:2)')';
compass=-atan2(csy(:,2),csy(:,1));
yaw = fmod(pi/2 - compass + 3*pi ,2*pi) - pi;

ts = (time - time(1))*24*3600;
I=find((ts>40) & (ts < 550));

figure(1);plot3(M(I,1),M(I,2),M(I,3),'+');
figure(2);plot3(Mr(I,1),Mr(I,2),Mr(I,3),'+');

Bned = B; Benu = [B(2), B(1), -B(3)];
Mr_ned = Mr;
Mr_enu = [Mr(:,2) Mr(:,1) -Mr(:,3)];

theta = linspace(0,2*pi,72);
idx = [100, 200, 400, 530]
for i=1:4;
    [v,k] = min(abs(ts - idx(i)))
    e = zeros(size(theta));
    for t=1:length(theta)
        e(t) = Mr_enu(k,:) * rotz(pi/2)*rotz(theta(t))' * Benu';
    end
    figure(i); plot(theta,e);
end
figure(5);plot(ts,yaw);


