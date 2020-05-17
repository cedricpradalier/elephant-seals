clear all

indivlist=['14686';'14876';'14878'; '14875';'14881';'14905';'15051';'15061'; '14873'; '14943'; '14676'; '14781'; '14899'; '14900'; '14901'; '14902'; '14903'; '14904'];

for i=1:size(indivlist,1)
	close all
	individual=indivlist(i,:)
	% dives=load('ml18_294dsens5.dives');
	% seg=load('ml18_294dsens5.seg');
	dir=['output.' individual]; 
	dives=load([dir '/divesum.csv']);
	seg=load([dir '/seg.csv']);
	dives(:,2) = dives(:,2);
	% dives=load('ml17_301asens5.dives');
	% seg=load('ml17_301asens5.seg');


	N=size(dives,1);
	% idx=[1:N]';
	idx=find((dives(:,2)>=1).*(dives(:,2)<=9));
	t=dives(idx,1)/86400;
	a=dives(idx,2);
	P=dives(idx,3);
	lt=linspace(min(t),max(t),60)';
	la=linspace(0,10,100)';


	idx0=find(dives(:,2)<0.1);
	t0=dives(idx0,1)/86400;
	P0=dives(idx0,3);
	A1=[ones(size(idx0)) t0 t0.*t0];
	B=P0;
	D1=A1\B;
	err = A1*D1-B;
	disp('A0: quadratic RMSE')
	sqrt(sum(err.^2))
	1 - mean(err.^2)/var(P0)

	A2=[ones(size(idx0)) t0 ];
	B=P0;
	D2=A2\B;
	err = A2*D2-B;
	disp('A0: linear RMSE')
	sqrt(sum(err.^2))
	1 - mean(err.^2)/var(P0)

	segd=seg(find((seg(:,4)==5).*((seg(:,8)-seg(:,6))>180.)),:);
	ts=segd(:,6)/86400;
	vs=(segd(:,9)-segd(:,7))./(segd(:,8)-segd(:,6));
	clf
	% plot(dives(idx0,1)/86400,dives(idx0,3),'+b',lt,[ones(size(lt)) lt lt.*lt]*D1,'-+r')
	subplot(2,1,1)
	plot(dives(idx0,1)/86400,dives(idx0,3),'.b',lt,[ones(size(lt)) lt lt.*lt]*D1,'-+r',lt,[ones(size(lt)) lt ]*D2,'-+g');
	grid on
	axis([0 60 -3 0])
	xlabel('t (j)')
	ylabel('max dP/dt (m/s)')
	title(sprintf('Individual %s Dive model 2D %.3f m/s/j',individual,D2(2)))
	subplot(2,1,2)
	plot(ts,vs,'*r');
	grid on
	axis([0 60 -0.6 0])
	xlabel('t (j)')
	ylabel('drift dP/dt (m/s)')
	print([dir '/divemodel2d.png'],'-dpng')



	A=[ones(size(idx)) t];
	B=P;
	C1=A\B;
	err = A*C1-B;
	disp('time only RMSE')
	sqrt(sum(err.^2))
	1 - mean(err.^2)/var(P)

	A=[ones(size(idx)) t a];
	B=P;
	C2=A\B;
	err = A*C2-B;
	disp('linear RMSE')
	sqrt(sum(err.^2))
	1 - mean(err.^2)/var(P)

	A=[ones(size(idx)) t a t.*t];
	B=P;
	C3=A\B;
	err = A*C3-B;
	disp('a-linear t-quadratic RMSE')
	sqrt(sum(err.^2))
	1 - mean(err.^2)/var(P)

	A=[ones(size(idx)) t a a.*a];
	B=P;
	C4=A\B;
	err = A*C4-B;
	disp('t-linear a-quadratic RMSE')
	sqrt(sum(err.^2))
	1 - mean(err.^2)/var(P)

	A=[ones(size(idx)) t a a.^a t.*t];
	B=P;
	C5=A\B;
	err = A*C5-B;
	disp('independent quadratic RMSE')
	sqrt(sum(err.^2))
	1 - mean(err.^2)/var(P)

	A=[ones(size(idx)) t a a.^a t.*t a.*t];
	B=P;
	C6=A\B;
	err = A*C6-B;
	disp('quadratic RMSE')
	sqrt(sum(err.^2))
	1 - mean(err.^2)/var(P)

	[Gt,Ga]=meshgrid(lt,la);
	gt=reshape(Gt,prod(size(Gt)),1);
	ga=reshape(Ga,prod(size(Ga)),1);
	% Am=[ones(size(gt)) gt ga ga.*ga gt.*gt ga.*gt];
	% Am=[ones(size(gt)) gt ga ga.*ga];
	Am=[ones(size(gt)) gt ga ];
	gdP=Am*C2;
	gdP=reshape(gdP,length(la),length(lt));
	clf
	plot3(dives(:,1)/86400,dives(:,2),dives(:,3),'.b');
	grid on
	hold on; mesh(Gt,Ga,gdP); hold off
	xlabel('t (j)')
	ylabel('az (m/s2)')
	zlabel('dP/dt (m/s)')
	title(sprintf('Individual %s Dive model 3D %.3f m/s/j %.3f m/s/az',individual,C2(2),C2(3)))
	print([dir '/divemodel3d.png'],'-dpng')

	% pkg load image
	load([dir '/spectro.dat']);
	maxP=spectro(1,2:end);
	t=spectro(2,2:end);
	f=spectro(2:end,1);
	S=spectro(3:end,2:end);
	% Ss=imsmooth(S,'Gaussian',5);
	Ss=imgaussfilt(S,5);
	% imagesc(t/86400,f,Ss)
	T=1./f;
	ifreq=find((T > 1.0).*(T < 3.0));
	Ssf=Ss(ifreq,:);

	F=f(ifreq)*ones(size(t));
	mean_f = sum(F .* Ssf)./sum(Ssf);
	mean_f2 = sum(F .* F .* Ssf)./sum(Ssf);
	var_f = mean_f2 - mean_f.*mean_f;
	std_f = sqrt(var_f);



	% Spectrogram over selected frequencies
	clf
	imagesc(t/86400,f(ifreq),Ssf)
	xlabel('t (j)')
	ylabel('f (Hz)')
	title(sprintf('Individual %s Climb spectrogram',individual))
	print([dir '/spectro.png'],'-dpng')
	% Correlation between frequency power and max depth
	clf
	plot(-maxP,max(Ssf),'+')
	ylabel('fft power')
	xlabel('-maxDepth')
	grid on
	title(sprintf('Individual %s Climb depth to fft power correlation',individual))
	print([dir '/depth2power.png'],'-dpng')

	% Correlation between swimming frequencies and maxP
	clf
	plot(-maxP,1./mean_f,'+')
	ylabel('swim period (s)')
	xlabel('-maxDepth')
	grid on
	title(sprintf('Individual %s Climb depth to swim period correlation',individual))
	print([dir '/depth2period.png'],'-dpng')


	% Co-evolution of frequency power and max depth
	clf
	plot(t/86400,max(Ssf),'-b')
	set(gcf,'PaperUnits','inches','PaperPosition',[0 0 60 10])
	ylabel('max fft power')
	xlabel('t (j)')
	grid on
	title(sprintf('Individual %s fft power evolution',individual))
	print([dir '/t2power.png'],'-dpng','-r100')

	clf
	plot(t/86400,std_f,'-r')
	set(gcf,'PaperUnits','inches','PaperPosition',[0 0 60 10])
	ylabel('std(f) (Hz)')
	xlabel('t (j)')
	grid on
	title(sprintf('Individual %s std(f) evolution',individual))
	print([dir '/t2stdf.png'],'-dpng','-r100')

	clf
	plot(t/86400,-maxP,'-k')
	set(gcf,'PaperUnits','inches','PaperPosition',[0 0 60 10])
	ylabel('-maxP (m)')
	xlabel('t (j)')
	grid on
	title(sprintf('Individual %s dive depth evolution',individual))
	print([dir '/t2maxP.png'],'-dpng','-r100')
end
