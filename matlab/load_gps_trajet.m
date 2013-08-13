fid = fopen('Trajet.csv', 'rt');
a = textscan(fid, '%d %s %f %f %f %f %f %f %f %f', 'Delimiter',',', 'CollectOutput',1, 'HeaderLines',1);
fclose(fid);

format short g
% M = [double(a{1}) datenum(a{2}) a{3}];
gps = [datenum(a{2}) a{3}(:,1:2) zeros(size(a{3},2))];
disp 'loaded GPS'

