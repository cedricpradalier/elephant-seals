clear all

% filelist=[ '2017-26-14622'; '2017-29-14331'; '2018-32-14875'; '2018-33-14873'; '2018-34-14879'; '2018-37-14876'; '2018-38-14881'; '2018-40-14878'; '2018-41-14686';...
% '2018-43-14781'; '2018-44-14904'; '2018-45-14900'; '2018-46-14901'; '2018-47-14905'; '2018-48-14676'; '2018-49-14903'; '2018-50-14899'; '2018-51-14902';...
% '2019-02-14899'; '2019-03-14902'; '2019-04-14903'; '2019-05-14901'; '2019-14-14873'; '2019-15-14876'; '2019-17-14875'; '2019-19-15061'; '2019-20-14878';...
% '2019-21-14881'; '2019-23-14905'; '2019-24-15051'];
filelist=[  '2019-12-14943' ; '2018-39-14874' ];
%filelist=['14686';'14876';'14878'; '14875';'14881';'14905';'15051';'15061'; '14873'; '14943'; '14676'; '14781'; '14899'; '14900'; '14901'; '14902'; '14903'; '14904'];
for i=1:size(filelist,1)
	filename=filelist(i,:)

	load([filename '_accelero.mat']);
	[t, index] = unique(acc_ac(:,1));
	A=acc_ac(index,2:4)*10; % Scaling in m/s2
	clear acc_ac index


	[tp,index]=unique(tdrcor2(:,1));
	P=tdrcor2(index,2);
	P=interp1(tp,P,t);
	clear index;

	%conversion in seconds
	t = (t-t(1)) * 86400;
	tr = [0:.2:max(t)]'; 
	A=[tr interp1(t,[A P],tr)];
	% A=resample([t A P],t,5);
	% save([filename '_tAP.mat'],'A');

	Aint=int32(A*100);
	disp 'finished loading'
	fid = fopen([filename '.csv'],'wt');
	fprintf(fid,'%d %d %d %d %d\n',Aint');
	fclose(fid);
	clear Aint;
	disp 'finished conversion to CSV'

	t=A(:,1);
	P=-A(:,5);
	A=A(:,2:4);

	HUNT=1;
	SURFACE=2;
	DIVE=3;
	CLIMB=4;
	DRIFT=5;

	label=ones(size(P)) * HUNT;
	% An=sqrt(sum((A.*A)')');
	% Surface
	idxS=find(P>-10);
	label(idxS)=SURFACE;

	Pf=imfilter(P,fspecial('average', [50 1]),'replicate');
	dPf=[0;diff(Pf)/0.2];
	idxS=find((label~=SURFACE).*(dPf<-1.0));
	label(idxS)=DIVE;
	idxS=find((label~=SURFACE).*(dPf>+0.5));
	label(idxS)=CLIMB;
	idxS=find((label~=SURFACE).*(dPf<-0.2).*(dPf>=-1.0).*(A(:,3)<-5));
	label(idxS)=DRIFT;
	disp 'completed first labeling'

	label=int8(label);
	%dlmwrite([filename '.labels'],label)
	fid = fopen([filename '.labels'],'wt');
	fprintf(fid,'%d\n',label');
	fclose(fid);

	disp 'Saved intermediate data'

end
