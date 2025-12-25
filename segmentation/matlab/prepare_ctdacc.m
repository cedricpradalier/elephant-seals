clear all
filelist=[ '14686'; '14873'; '14875'; '14876'; '14878';...
'14881'; '14905'; '14943'; '15051'; '15061']

for i=1:size(filelist,1)
	filename=filelist(i,:)
	load([filename '_tAP.mat']);
	Aint=int32(A*100);
	disp 'finished loading'
	fid = fopen([filename '.csv'],'wt');
	fprintf(fid,'%d %d %d %d %d\n',Aint');
	fclose(fid);

	%csvwrite([filename '.csv'],Aint);
	%save([filename '.csv'],'Aint','-ascii');
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
	%dlmwrite([filename '.labels'],label)
	fid = fopen([filename '.labels'],'wt');
	fprintf(fid,'%d\n',label');
	fclose(fid);

	disp 'Saved intermediate data'
end
return
