clear all
% close all

addpath('rotations');
addpath('optimisation');

load('../horses.txt');
data = horses;
mag_offset = 0.5;

N = size(data,1)
A = data(:,2:4);
M = data(:,5:7);
nA=sqrt(sum((A.*A)')');
nM=sqrt(sum((M.*M)')');
An = A .* ((1./nA) * [1 1 1]);
Mn = M .* ((1./nM) * [1 1 1]);
ts = data(:,1);

% - sign due to gravity being negative
pitch=asin(An(:,1));
roll=-atan2(An(:,2),-An(:,3));
% Mr = zeros(N,3);
% for i=1:N
%     Mr(i,:) = (rpy(roll(i),pitch(i),0) * Mn(i,:)')';
% end
% csy=(Bmat * Mr(:,1:2)')';
% compass=-atan2(csy(:,2),csy(:,1));
yaw = atan2(Mn(:,2),Mn(:,1)) - mag_offset;

figure(1);plot(ts,roll*180/pi,'b-');grid on
print -dpng 'horse_roll.png'
figure(2);plot(ts,pitch*180/pi,'r-');grid on
print -dpng 'horse_pitch.png'
figure(3);plot(ts,yaw*180/pi,'g+');grid on
print -dpng 'horse_yaw.png'

figure(4);plot(ts,An);grid on
print -dpng 'horse_acc.png'
figure(5);plot(ts,Mn);grid on
print -dpng 'horse_mag.png'

