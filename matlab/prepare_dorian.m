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


load('data/Data_acousonde626019.mat')

disp 'loaded core data'

M = Magneto_acousonde626019;
A = Accelero_acousonde626019;
N = size(A,1)
nA=sqrt(sum((A.*A)')');
nM=sqrt(sum((M.*M)')');
An = A .* ((1./nA) * [1 1 1]);
Mn = M .* ((1./nM) * [1 1 1]);


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

D=ones(N,1);

preload = [cumsum(D), cumsum(D), A, M, D, D, D, D,D,D, roll, pitch, yaw, D, D];
save -ascii data/preload_mat_626019.txt preload

load('data/Data_acousonde626040.mat')


disp 'loaded core data'

M = Magneto_acousonde626040;
A = Accelero_acousonde626040;
N = size(A,1)
nA=sqrt(sum((A.*A)')');
nM=sqrt(sum((M.*M)')');
An = A .* ((1./nA) * [1 1 1]);
Mn = M .* ((1./nM) * [1 1 1]);


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

D=ones(N,1);

preload = [cumsum(D), cumsum(D), A, M, D, D, D, D,D,D, roll, pitch, yaw, D, D];
save -ascii data/preload_mat_626040.txt preload

