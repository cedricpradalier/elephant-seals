clear all
close all

addpath('rotations');
addpath('optimisation');

load('data/BioVar626019_5Hz.mat')
idx=find(~isnan(Latitude));
Latitude=Latitude(idx);
Longitude=Longitude(idx);
TimeVector=TimeVector(idx);
clear Accelero_x Accelero_y Accelero_z
clear Magneto_x Magneto_y Magneto_z
clear Depth StatusDive Roll Pitch Yaw
B = []
addpath('GPS_CoordinateXforms')
addpath('igrf')

ngps=size(Latitude,1);
[x,y,utmzone] = wgs2utm(Latitude,Longitude,43);
B = ones(ngps,1) * igrf(TimeVector(1),Latitude(1),Longitude(1),0);
% for i=1:1000000:ngps
%     B = [B;igrf(TimeVector(1),Latitude(i:min(i+1000000-1,ngps)),...
%         Longitude(i:min(i+1000000-1,ngps)),0)];
% end
% Outputs:
%   -BX: Northward component of the magnetic field in nanoteslas (nT).
%   -BY: Eastward component of the magnetic field in nT.
%   -BZ: Downward component of the magnetic field in nT.
%   -B: [BX(:), BY(:), BZ(:)].
B = B ./ (sqrt(sum((B.*B)')')*ones(1,3));
disp 'Got magnetic field data'
gps2 = [floor(TimeVector(:,1)) (TimeVector(:,1)-floor(TimeVector(:,1))) ...
    Longitude Latitude Longitude Latitude x y utmzone B];
save -ascii data/gps_mat_all.txt gps2
Bmat = inv([B(1,1) -B(1,2);B(1,2) B(1,1)]);
clear B x y utmzone Longitude Latitude gps2

disp 'loaded core data'
load('data/BioVar626019_5Hz.mat')
clear Longitude Latitude StatusDive Depth

M = [Magneto_x(idx) Magneto_y(idx) Magneto_z(idx)];
A = [Accelero_x(idx) Accelero_y(idx) Accelero_z(idx)];
M = M(1:100000,:);
A = A(1:100000,:);
N = size(M,1)
clear Accelero_x Accelero_y Accelero_z
clear Magneto_x Magneto_y Magneto_z
if true
    % clear Roll Pitch Yaw
    nA=sqrt(sum((A.*A)')');
    nM=sqrt(sum((M.*M)')');
    An = A .* ((1./nA) * [1 1 1]);
    Mn = M .* ((1./nM) * [1 1 1]);


    % - sign due to gravity being negative
    pitch=asin(An(:,1));
    roll=-atan2(An(:,2),-An(:,3));
    % Ar = zeros(N,3);
    % target=0.05;
    % for i=1:N
    %     if i/N > target
    %         target = target + 0.05
    %         disp(sprintf('Ar: %d%% done',i/N))
    %     end
    %     Ar(i,:) = (rpy(roll(i),pitch(i),0) * A(i,:)')';
    % end
    % disp('Done Ar')
    target=0.05;
    Mr = zeros(N,3);
    for i=1:N
        if i/N > target
            target = target + 0.05
            disp sprintf('Mr: %d%% done',i/N)
        end
        Mr(i,:) = (rpy(roll(i),pitch(i),0) * Mn(i,:)')';
    end
    disp('Done Mr')
    csy=(Bmat * Mr(:,1:2)')';
    cps=-atan2(csy(:,2),csy(:,1));
    yaw = -mod(pi/2 - cps ,2*pi) ;

    disp 'Extracted RPY'
end

D=ones(N,1);

preload = [floor(TimeVector(idx,1)) (TimeVector(idx,1)-floor(TimeVector(idx,1))) ...
    A, M, D, D, D, D,D,D, Roll(idx), Pitch(idx), Yaw(idx), D, D];
save -ascii data/BioVar626019_5Hz.txt preload

