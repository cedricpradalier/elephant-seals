clear all
close all

addpath('rotations');
addpath('optimisation');

load('data/gps.mat');
if size(gps,2)==5
    B = []
    addpath('GPS_CoordinateXforms')
    addpath('igrf')

    ngps=size(gps,1);
    [x,y,utmzone] = wgs2utm(gps(:,3),gps(:,2),43);
    for i=1:ngps
        B = [B;igrf(data(1,1),gps(i,3),gps(i,2),0)];
    end
    % Outputs:
    %   -BX: Northward component of the magnetic field in nanoteslas (nT).
    %   -BY: Eastward component of the magnetic field in nT.
    %   -BZ: Downward component of the magnetic field in nT.
    %   -B: [BX(:), BY(:), BZ(:)].
    B = B ./ (sqrt(sum((B.*B)')')*ones(1,3));
    disp 'Got magnetic field data'
    gps = [gps x y utmzone B]
end
B = gps(:,9:11);
Bmat = inv([B(1,1) -B(1,2);B(1,2) B(1,1)]);


load('data/dataspd.mat');
% gps = gps(:,1:3);
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


% Select a single cycle
idx = find(gps(:,4)>0);
data = acou(gps(idx(1),4):gps(idx(2),4),:);
gps = gps(idx(1):idx(2),:);
gps(:,4) = gps(:,4) - gps(1,4) + 1;
disp 'loaded core data'

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

dive = data(:,13);
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
data(:,13) = dive;
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
C=[sum(W.*x.*x) sum(W.*x);sum(W.*x) length(x)];
B=[sum(W.*x.*y);sum(W.*y)];
p=inv(C)*B;
axfit=sin(p(1)*qpbin + p(2));
% scatter(qpitch,An(:,1),2,dive);caxis([-0.5 4.5]); colorbar; axis equal; grid on;
% hold on;
% plot(qpbin,axmean,'-r',qpbin,axfit,'-g','linewidth',4); 
% hold off
disp 'Prepared speed model'



% - sign due to gravity being negative
pitch=asin(An(:,1));
roll=-atan2(An(:,2),-An(:,3));
Ar = zeros(N,3);
for i=1:N
    Ar(i,:) = (rpy(roll(i),pitch(i),0) * A(i,:)')';
end
Mr = zeros(N,3);
% Anone = zeros(N,3);
%Anone2 = zeros(N,3);
for i=1:N
    Mr(i,:) = (rpy(roll(i),pitch(i),0) * Mn(i,:)')';
end
csy=(Bmat * Mr(:,1:2)')';
compass=-atan2(csy(:,2),csy(:,1));
yaw = -mod(pi/2 - compass ,2*pi) ;

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
for i=1:N-1
    R = rpy(roll(i),pitch(i),yaw(i));
    vdt(i+1,:) = R * [ds(i), 0, 0]';
end

X=cumsum(vdt);X(:,3) = depth;

day = floor(data(:,1));
preload = [day, (data(:,1)-day), data(:,7:9), data(:,4:6), depth, vel, has_vel, X, roll, pitch, yaw, vel_all, dive];
preload_titles ={};
preload_titles{1} = 'datenum';
preload_titles{2} = 'datenum2';
preload_titles{3} = 'Ax';
preload_titles{4} = 'Ay';
preload_titles{5} = 'Az';
preload_titles{6} = 'Mx';
preload_titles{7} = 'My';
preload_titles{8} = 'Mz';
preload_titles{9} = 'depth';
preload_titles{10} = 'vel';
preload_titles{11} = 'has_vel';
preload_titles{12} = 'X';
preload_titles{13} = 'Y';
preload_titles{14} = 'Z';
preload_titles{15} = 'roll';
preload_titles{16} = 'pitch';
preload_titles{17} = 'yaw';
preload_titles{18} = 'V';   
preload_titles{19} = 'dive_status';   
save -mat data/preload_trajet.mat preload preload_titles
save -ascii data/preload_mat.txt preload
gps2 = [floor(gps(:,1)) (gps(:,1)-floor(gps(:,1))) gps(:,2:end)];
save -ascii data/gps_mat.txt gps2

Mb = zeros(size(Mn));
for i=1:N
    R = rpy(roll(i),pitch(i),-yaw(i));
    Mb(i,:) = (R * Mn(i,:)')';
end


