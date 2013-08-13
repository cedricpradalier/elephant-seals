clear all
close all
addpath('../..');
fid = fopen('Trajet.csv', 'rt');
a = textscan(fid, '%d %s %f %f %f %f %f %f %f %f', 'Delimiter',',', 'CollectOutput',1, 'HeaderLines',1);
fclose(fid);

format short g
% M = [double(a{1}) datenum(a{2}) a{3}];
gps = [datenum(a{2}) a{3}(:,1:2) zeros(size(a{3},2))];
disp 'loaded GPS'

load('dataspd.mat');

gps = gps(:,1:3);
acou_surface=acou(find(acou(:,13)==0),:);
for i=1:size(gps,1)
    [v,w] = min(abs(acou_surface(:,1)-gps(i,1)));
    w = find(acou(:,1)==acou_surface(w,1));
    gps(i,5) = v;
    if abs(v)<(30./(24*60))
        gps(i,4) = w;
    else
        gps(i,4) = 0;
    end
end

% P =[   0.717812 ;    0.050772 ]; % For asin model
P =[ 0.9242688148141405; 0.0839373429529526]; % for linear model
P =[1.113828;  -0.011317];
data = acou;
disp 'loaded core data'

B = []
addpath('../../GPS_CoordinateXforms')
addpath('../../igrf')

ngps=size(gps,1);
for i=1:ngps
    B = [B;igrf(data(1,1),gps(i,3),gps(i,2),0)];
end
% Outputs:
%   -BX: Northward component of the magnetic field in nanoteslas (nT).
%   -BY: Eastward component of the magnetic field in nT.
%   -BZ: Downward component of the magnetic field in nT.
%   -B: [BX(:), BY(:), BZ(:)].
B = B ./ (sqrt(sum((B.*B)')')*ones(1,3));
Bmat = inv([B(1,1) -B(1,2);B(1,2) B(1,1)]);
disp 'Got magnetic field data'

N = size(data,1)
vel=data(:,22);
has_vel=data(:,21);
M = data(:,4:6);
A = data(:,7:9);
nA=sqrt(sum((A.*A)')');
nM=sqrt(sum((M.*M)')');
An = A .* ((1./nA) * [1 1 1]);
Mn = M .* ((1./nM) * [1 1 1]);
ts = (data(:,1) - data(1,1)) * 3600 * 24;

% dive=-ones(N,1);
% dive(find(data(:,1)>-5)) = 0;
% dive(find(data(:,1)<-300)) = 2;
% dive(40821)=-1;
% ndive=[dive(2:end);-1];
% takingoff=find((ndive==-1)&(dive==2));
% landing=find((ndive==2)&(dive==-1));
% surfacing=find((ndive==0)&(dive==-1));
% diving=find((ndive==-1)&(dive==0));
% for i=1:10
%     dive(takingoff(i):surfacing(i))=3;
%     dive(diving(i):landing(i))=1;
% end
dive = data(:,13);
startdive=find((dive(1:N-1)==0)&(dive(2:N)==1));
depth = data(:,2);
ddepth=[0;diff(depth)];
dt=[0;diff(ts)];
disp 'Affected default columns and normalised'


dzdt=ddepth./dt;
velf = vel;
velf(find(dive==0))=0;
qpitch=(dzdt./velf);
qpitch(find((qpitch)>1))=1;
qpitch(find((qpitch)<-1))=-1;
vpitch=-asin(qpitch);

dive(find((dive~=0)&(dive~=3)&(qpitch<-0.5)&(An(:,1)>0))) = 4;
disp 'prepared dive coefficients'

delta=0.05;
qpbin=-1:delta:1;
J = find(abs(qpitch)~=1);
qp = qpitch(J);
ap = An(J,1);
axmean = zeros(size(qpbin));
W = zeros(size(qpbin));
for i=1:length(qpbin)
    I = find((abs(qp-qpbin(i))<delta/2)&(dive(J)~=4));
    W(i) = length(I)/N;
    axmean(i)=mean(ap(I));
end
% plot(qpbin,axmean); grid on;
x=qpbin;y=asin(axmean);
A=[sum(W.*x.*x) sum(W.*x);sum(W.*x) length(x)];
B=[sum(W.*x.*y);sum(W.*y)];
p=inv(A)*B;
axfit=sin(p(1)*qpbin + p(2));
% scatter(qpitch,An(:,1),2,dive);caxis([-0.5 4.5]); colorbar; axis equal; grid on;
% hold on;
% plot(qpbin,axmean,'-r',qpbin,axfit,'-g','linewidth',4); 
% hold off
disp 'Prepared speed model'



% - sign due to gravity being negative
pitch=-asin(An(:,1));
% Aroll = zeros(N,3);
% for i=1:N
%     R = rpy(0,-pitch(i),0);
%     Aroll(i,:) = (R * An(i,:)')';
% end
% roll=asin(Aroll(:,2));
roll=atan2(An(:,2),-An(:,3));
Mr = zeros(N,3);
% Anone = zeros(N,3);
Anone2 = zeros(N,3);
for i=1:N
    %R = rpy(-roll(i),0,0);
    %Anone(i,:) = (R * Aroll(i,:)')';
    R = rpy(0,-pitch(i),0) * rpy(-roll(i),0,0);
    Anone2(i,:) = (R * An(i,:)')';
    Mr(i,:) = (R * Mn(i,:)')';
end
csy=(Bmat * Mr(:,1:2)')';
compass=-atan2(csy(:,2),csy(:,1));
yaw = mod(pi/2 - compass ,2*pi) ;

disp 'Extracted RPY'
% Mnone = zeros(N,3);
% for i=1:N
%     R = rpy(0,0,-yaw(i));
%     Mnone(i,:) = (R * Mr(i,:)')';
% end
% yaw=pi - yaw
% plot(Mnone(:,1),Mnone(:,2),'.')


% vpred=-dzdt ./ sin(P(1)*pitch+P(2)); % linear model
% vpred=-dzdt ./ (P(1)*pitch+P(2)); % asin model
vpred = dzdt*p(1) ./ (asin(An(:,1)) - p(2));
vpred=min(3,max(0,vpred));
wndw=10;
vpredf=conv(vpred,ones(wndw,1)/wndw);vpredf=vpredf((wndw/2):(N+wndw/2-1));
vpredf(find(dive==0)) = 0;

vel_all = vel.*has_vel + vpredf.*(1-has_vel);
disp 'Computed velocity prediction'

ds=vel_all.*dt.*(dive~=0);

vdt=zeros(N,3);
for i=2:N
    R = rpy(roll(i),pitch(i),yaw(i));
    vdt(i,:) = R * [ds(i), 0, 0]';
end

X=cumsum(vdt);X(:,3) = depth;

J=find(dive==0);
igps = zeros(ngps,3);
for i=1:ngps
    [val,idx] = min(abs(data(J,1)-gps(i,1)));
    igps(i,1) = J(idx);
end
Y=[X(:,1)-X(igps(1,1),1)+x(1) X(:,2)-X(igps(1,1),2)+y(1) X(:,3)];
[x,y,utmzone] = wgs2utm(gps(:,3),gps(:,2),43);
igps(:,2:3) = [x - Y(igps(:,1),1), y - Y(igps(:,1),2)];
disp 'Computed integrated trajectory'
plot(Y(:,1),Y(:,2)-4e6,'-',Y(J,1),Y(J,2)-4e6,'*r',x,y-4e6,'og','markersize',20)
hold on
quiver(Y(igps(:,1),1),Y(igps(:,1),2)-4e6,igps(:,2),igps(:,3),0,'m')
hold off
axis equal
  
plot3(Y(:,1),Y(:,2)-4e6,Y(:,3),'-',Y(:,1),Y(:,2)-4e6,-550*ones(N,1),'-b',Y(J,1),Y(J,2)-4e6,-550*ones(size(J)),'*r',x,y-4e6,-550*ones(4,1),'og','markersize',10)
% axis([917500 919500 516000 519000 -560 10])


% I=1:10:N;
% I=find((dive~=0)&(dive~=2));
% figure(1); scatter(A(I,1),A(I,2),4,dive(I));caxis([-0.5 3.5]); colorbar; axis equal; grid on;
% figure(2); scatter(A(I,2),A(I,3),4,dive(I));caxis([-0.5 3.5]); colorbar; axis equal; grid on;
% figure(3); scatter(A(I,1),A(I,3),4,dive(I));caxis([-0.5 3.5]); colorbar; axis equal; grid on;
% 
% figure(1); scatter(M(I,1),M(I,2),4,dive(I));colorbar; axis equal;
% figure(2); scatter(M(I,2),M(I,3),4,dive(I));colorbar; axis equal;
% figure(3); scatter(M(I,1),M(I,3),4,dive(I));colorbar; axis equal;
% 
% 
% figure(1); scatter(Mr(I,1),Mr(I,2),4,dive(I));caxis([-0.5 3.5]); colorbar; axis equal; grid on;
% figure(2); scatter(Mr(I,2),Mr(I,3),4,dive(I));caxis([-0.5 3.5]); colorbar; axis equal; grid on;
% figure(3); scatter(Mr(I,1),Mr(I,3),4,dive(I));caxis([-0.5 3.5]); colorbar; axis equal; grid on;


