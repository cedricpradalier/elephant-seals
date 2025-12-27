clear all

indivlist={}
indivlist{1,1}=[ '2017-26-14622'; '2017-29-14331'; '2018-32-14875'; '2018-33-14873'; '2018-34-14879'; '2018-37-14876'; '2018-38-14881'; '2018-40-14878'; '2018-41-14686';...
'2018-43-14781'; '2018-44-14904'; '2018-45-14900'; '2018-46-14901'; '2018-47-14905'; '2018-48-14676'; '2018-49-14903'; '2018-50-14899'; '2018-51-14902';...
'2019-02-14899'; '2019-03-14902'; '2019-04-14903'; '2019-05-14901'; '2019-14-14873'; '2019-15-14876'; '2019-17-14875'; '2019-19-15061'; '2019-20-14878';...
'2019-21-14881'; '2019-23-14905'; '2019-24-15051'; '2019-12-14943' ; '2018-39-14874'];
indivlist{1,2}=['ml17301';'ml18294'];
indivlist{1,3}=['15734-inter'];

% for ic=1:size(indivlist,2)
for ic=3:3
	for i=1:size(indivlist{1,ic},1)
		close all
		individual=indivlist{1,ic}(i,:)
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
		idx=find((dives(:,2)>=1).*(dives(:,2)<=9).*(dives(:,3)>-5));
		t=dives(idx,1)/86400;
		a=dives(idx,2);
		P=dives(idx,3);
		lt=linspace(min(t),max(t),60)';
		la=linspace(0,10,100)';


		idx0=find((dives(:,2)<0.1).*(dives(:,3)>-5));
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

		segd=seg(find((seg(:,5)==5).*((seg(:,9)-seg(:,7))>180.)),:);
		ts=segd(:,7)/86400;
		vs=(segd(:,10)-segd(:,8))./(segd(:,9)-segd(:,7));
		clf
		% plot(dives(idx0,1)/86400,dives(idx0,3),'+b',lt,[ones(size(lt)) lt lt.*lt]*D1,'-+r')
		subplot(2,1,1)
		plot(dives(idx0,1)/86400,dives(idx0,3),'.b',lt,[ones(size(lt)) lt lt.*lt]*D1,'-+r',lt,[ones(size(lt)) lt ]*D2,'-+g');
		grid on
		A=axis()
		axis([0 A(2) -3 0])
		xlabel('t (j)')
		ylabel('max dP/dt (m/s)')
		title(sprintf('Individual %s Dive model 2D %.3f m/s/j',individual,D2(2)))
		subplot(2,1,2)
		plot(ts,vs,'*r');
		grid on
		axis([0 A(2) -0.6 0])
		xlabel('t (j)')
		ylabel('drift dP/dt (m/s)')
		print([dir '/divemodel2d.png'],'-dpng')

		driftrate=[ts*86400 vs];
		save([dir '/driftrate.csv'],'-ascii','-double','driftrate');


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
		A=axis();
		axis([A(1) A(2) 0 10 -3 0]);
		title(sprintf('Individual %s Dive model 3D %.3f m/s/j %.3f m/s/az',individual,C2(2),C2(3)))
		print([dir '/divemodel3d.png'],'-dpng')

		close all
		% pkg load image
		load([dir '/spectro_climb.dat']);
		spectro=spectro_climb;
		clear spectro_climb;
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
		print([dir '/climb_spectro.png'],'-dpng')
		% Correlation between frequency power and max depth
		clf
		plot(-maxP,max(Ssf),'+')
		ylabel('fft power')
		xlabel('-maxDepth')
		grid on
		title(sprintf('Individual %s Climb depth to fft power correlation',individual))
		print([dir '/climb_depth2power.png'],'-dpng')

		% Correlation between swimming frequencies and maxP
		clf
		plot(-maxP,1./mean_f,'+')
		ylabel('swim period (s)')
		xlabel('-maxDepth')
		grid on
		title(sprintf('Individual %s Climb depth to swim period correlation',individual))
		print([dir '/climb_depth2period.png'],'-dpng')

		A=[ones(size(t))' t'/86400 -maxP']; B=max(Ssf)';
		X=A\B;
		err = A*X-B;
		disp('Linear relation between max(Ssf) and (t,maxP)')
		sqrt(sum(err.^2))
		1 - mean(err.^2)/var(B)

		lt=linspace(min(t/86400),max(t/86400),20)';
		lP=linspace(0,max(-maxP),20);
		[Gt,GP]=meshgrid(lt,lP);
		gt=reshape(Gt,prod(size(Gt)),1);
		gP=reshape(GP,prod(size(GP)),1);
		Am=[ones(size(gt)) gt gP ];
		gdP=Am*X;
		gdP=reshape(gdP,length(lP),length(lt));
		clf
		plot3(t/86400,-maxP,max(Ssf),'.');
		grid on
		hold on; mesh(Gt,GP,gdP); hold off
		xlabel('t (j)')
		ylabel('-maxP (m)')
		zlabel('fft power')
		title(sprintf('Individual %s Climb depth to fft power correlation, over time',individual))
		print([dir '/climb_depthtime2power.png'],'-dpng')

		% Co-evolution of frequency power and max depth
		clf
		plot(t/86400,max(Ssf),'-k','linewidth',2)
		set(gcf,'PaperUnits','inches','PaperPosition',[0 0 60 10])
		ylabel('max fft power')
		xlabel('t (j)')
		grid on
		title(sprintf('Individual %s fft power evolution',individual))
		print([dir '/climb_t2power.png'],'-dpng','-r100')

		clf
		plot(t/86400,std_f,'-k','linewidth',2)
		set(gcf,'PaperUnits','inches','PaperPosition',[0 0 60 10])
		ylabel('std(f) (Hz)')
		xlabel('t (j)')
		grid on
		title(sprintf('Individual %s std(f) evolution',individual))
		print([dir '/climb_t2stdf.png'],'-dpng','-r100')

		clf
		plot(t/86400,-maxP,'-k')
		set(gcf,'PaperUnits','inches','PaperPosition',[0 0 60 10])
		ylabel('-maxP (m)')
		xlabel('t (j)')
		grid on
		title(sprintf('Individual %s dive depth evolution',individual))
		print([dir '/climb_t2maxP.png'],'-dpng','-r100')

		clf
		idx=find(-maxP>200);
		lt=linspace(0,max(t)/86400,20)';
		A=[ones(size(idx))' t(idx)'/86400]; B=max(Ssf(:,idx))./(-maxP(idx));
		X=A\B';
		lB=[ones(size(lt)) lt]*X;
		plot(t(idx)/86400,max(Ssf(:,idx))./(-maxP(idx)),'.k',lt,lB,'-g','linewidth',2)
		set(gcf,'PaperUnits','inches','PaperPosition',[0 0 60 10])
		ylabel('fft power/maxP')
		xlabel('t (j)')
		grid on
		title(sprintf('Individual %s normalised fft power evolution',individual))
		print([dir '/climb_t2normpower.png'],'-dpng','-r100')

		is_good=zeros(size(t));
		is_good(idx)=ones(size(idx));
		climbfreqstats=[t' is_good' -maxP' std_f' max(Ssf)'];
		save([dir '/climbfreqstats.csv'],'-ascii','-double','climbfreqstats');

		close all
		% pkg load image
		load([dir '/spectro_dive.dat']);
		spectro=spectro_dive;
		clear spectro_dive;
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
		title(sprintf('Individual %s Dive spectrogram',individual))
		print([dir '/dive_spectro.png'],'-dpng')
		% Correlation between frequency power and max depth
		clf
		plot(-maxP,max(Ssf),'+')
		ylabel('fft power')
		xlabel('-maxDepth')
		grid on
		title(sprintf('Individual %s Climb depth to fft power correlation',individual))
		print([dir '/dive_depth2power.png'],'-dpng')

		% Correlation between swimming frequencies and maxP
		clf
		plot(-maxP,1./mean_f,'+')
		ylabel('swim period (s)')
		xlabel('-maxDepth')
		grid on
		title(sprintf('Individual %s Climb depth to swim period correlation',individual))
		print([dir '/dive_depth2period.png'],'-dpng')

		A=[ones(size(t))' t'/86400 -maxP']; B=max(Ssf)';
		X=A\B;
		err = A*X-B;
		disp('Linear relation between max(Ssf) and (t,maxP)')
		sqrt(sum(err.^2))
		1 - mean(err.^2)/var(B)

		lt=linspace(min(t/86400),max(t/86400),20)';
		lP=linspace(0,max(-maxP),20);
		[Gt,GP]=meshgrid(lt,lP);
		gt=reshape(Gt,prod(size(Gt)),1);
		gP=reshape(GP,prod(size(GP)),1);
		Am=[ones(size(gt)) gt gP ];
		gdP=Am*X;
		gdP=reshape(gdP,length(lP),length(lt));
		clf
		plot3(t/86400,-maxP,max(Ssf),'.');
		grid on
		hold on; mesh(Gt,GP,gdP); hold off
		xlabel('t (j)')
		ylabel('-maxP (m)')
		zlabel('fft power')
		title(sprintf('Individual %s Climb depth to fft power correlation, over time',individual))
		print([dir '/dive_depthtime2power.png'],'-dpng')

		% Co-evolution of frequency power and max depth
		clf
		plot(t/86400,max(Ssf),'-k','linewidth',2)
		set(gcf,'PaperUnits','inches','PaperPosition',[0 0 60 10])
		ylabel('max fft power')
		xlabel('t (j)')
		grid on
		title(sprintf('Individual %s fft power evolution',individual))
		print([dir '/dive_t2power.png'],'-dpng','-r100')

		clf
		plot(t/86400,std_f,'-k','linewidth',2)
		set(gcf,'PaperUnits','inches','PaperPosition',[0 0 60 10])
		ylabel('std(f) (Hz)')
		xlabel('t (j)')
		grid on
		title(sprintf('Individual %s std(f) evolution',individual))
		print([dir '/dive_t2stdf.png'],'-dpng','-r100')

		clf
		plot(t/86400,-maxP,'-k')
		set(gcf,'PaperUnits','inches','PaperPosition',[0 0 60 10])
		ylabel('-maxP (m)')
		xlabel('t (j)')
		grid on
		title(sprintf('Individual %s dive depth evolution',individual))
		print([dir '/dive_t2maxP.png'],'-dpng','-r100')

		clf
		idx=find(-maxP>200);
		lt=linspace(0,max(t)/86400,20)';
		A=[ones(size(idx))' t(idx)'/86400]; B=max(Ssf(:,idx))./(-maxP(idx));
		X=A\B';
		lB=[ones(size(lt)) lt]*X;
		plot(t(idx)/86400,max(Ssf(:,idx))./(-maxP(idx)),'.k',lt,lB,'-g','linewidth',2)
		set(gcf,'PaperUnits','inches','PaperPosition',[0 0 60 10])
		ylabel('fft power/maxP')
		xlabel('t (j)')
		grid on
		title(sprintf('Individual %s normalised fft power evolution',individual))
		print([dir '/dive_t2normpower.png'],'-dpng','-r100')

		is_good=zeros(size(t));
		is_good(idx)=ones(size(idx));
		divefreqstats=[t' is_good' -maxP' std_f' max(Ssf)'];
		save([dir '/divefreqstats.csv'],'-ascii','-double','divefreqstats');

		%disp('press enter to continue')
		%pause
	end
end
