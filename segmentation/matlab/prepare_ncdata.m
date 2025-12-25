clear all
%pkg load image
%pkg load netcdf

filelist=['ml22_291atrk']; % 'ml17_301asens5';'ml18_294dsens5'];
for i=1:size(filelist,1)
	filename=filelist(i,:)
	% filename='ml17_280asens5_6hrs'
	A=ncread([filename '.nc'],'A');
	P=-ncread([filename '.nc'],'P');
	% J=ncread(filename,'J');
	% LL=ncread(filename,'LL');
	% T=ncread(filename,'T');
	% M=ncread(filename,'M');
	disp 'finished loading'

	HUNT=1;
	SURFACE=2;
	DIVE=3;
	CLIMB=4;
	DRIFT=5;

	label=ones(size(P)) * HUNT;
	t=cumsum(ones(size(P))*0.2);
	An=sqrt(sum((A.*A)')');
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
	Aint=int32([t A -P]*100);
	fid = fopen([filename '.csv'],'wt');
	fprintf(fid,'%d %d %d %d %d\n',Aint');
	fclose(fid);
	fid = fopen([filename '.labels'],'wt');
	fprintf(fid,'%d\n',label');
	fclose(fid);


	%save('-ascii',[filename '.csv'],'label')
	%csvwrite([filename '.csv'],Aggregate,'precision','%.2f')
	%csvwrite([filename '.labels'],label)
	clear Aint

	disp 'Saved intermediate data'
end
