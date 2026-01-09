clear all

indivlist={}
indivlist{1,1}=[ '2017-26-14622'; '2017-29-14331'; '2018-32-14875'; '2018-33-14873'; '2018-34-14879'; '2018-37-14876'; '2018-38-14881'; '2018-40-14878'; '2018-41-14686';...
'2018-43-14781'; '2018-44-14904'; '2018-45-14900'; '2018-46-14901'; '2018-47-14905'; '2018-48-14676'; '2018-49-14903'; '2018-50-14899'; '2018-51-14902';...
'2019-02-14899'; '2019-03-14902'; '2019-04-14903'; '2019-05-14901'; '2019-14-14873'; '2019-15-14876'; '2019-17-14875'; '2019-19-15061'; '2019-20-14878';...
'2019-21-14881'; '2019-23-14905'; '2019-24-15051'; '2019-12-14943' ; '2018-39-14874'];
indivlist{1,2}=['ml17301';'ml18294'];

individ=0;
for ic=1:size(indivlist,2)
	for i=1:size(indivlist{1,ic},1)
		individ=individ + 1;


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
		% Select non-vertical motion. Unused
		idx=find((dives(:,2)>=1).*(dives(:,2)<=9).*(dives(:,3)>-5));
		t=dives(idx,1)/86400;
		a=dives(idx,2);
		P=dives(idx,3);

		% Select az < 0.1 (i.e. seal is vertical) and seal at least 5m from the surface
		idx0=find((dives(:,2)<0.1).*(dives(:,3)>-5));
		t0=dives(idx0,1)/86400;
		P0=dives(idx0,3);
		[t1,ia1,ic1] = unique(t0);
		pique=[t1,P0(ia1)];

		

		% Select drift and segment length > 3 min
		segd=seg(find((seg(:,5)==5).*((seg(:,9)-seg(:,7))>180.)),:);
		ts=segd(:,7)/86400;
		vs=(segd(:,10)-segd(:,8))./(segd(:,9)-segd(:,7));
		drift=[ts vs];

		% Can be used to provide a dive id to any measurement with interp1d
		time_and_dive=[seg(:,1) seg(:,7)/86400];


		% clf
		% plot(dives(idx0,1)/86400,dives(idx0,3),'+b',lt,[ones(size(lt)) lt lt.*lt]*D1,'-+r')
		% subplot(2,1,1)
		% plot(dives(idx0,1)/86400,dives(idx0,3),'.b');
		% grid on
		% A=axis()
		% axis([0 A(2) -3 0])
		% xlabel('t (j)')
		% ylabel('max dP/dt (m/s)')
		% title(sprintf('Individual %s Dive model 2D %.3f m/s/j',individual,D2(2)))
		% subplot(2,1,2)
		% plot(ts,vs,'*r');
		% grid on
		% axis([0 A(2) -0.6 0])
		% xlabel('t (j)')
		% ylabel('drift dP/dt (m/s)')
		% print([dir '/divemodel2d.png'],'-dpng')



		% clf
		% plot3(dives(:,1)/86400,dives(:,2),dives(:,3),'.b');
		% grid on
		% xlabel('t (j)')
		% ylabel('az (m/s2)')
		% zlabel('dP/dt (m/s)')
		% A=axis();
		% axis([A(1) A(2) 0 10 -3 0]);
		% title(sprintf('Individual %s Dive model 3D %.3f m/s/j %.3f m/s/az',individual,C2(2),C2(3)))
		% print([dir '/divemodel3d.png'],'-dpng')

		% close all
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
		% clf
		% imagesc(t/86400,f(ifreq),Ssf)
		% xlabel('t (j)')
		% ylabel('f (Hz)')
		% title(sprintf('Individual %s Climb spectrogram',individual))
		% print([dir '/climb_spectro.png'],'-dpng')
		% Correlation between frequency power and max depth
		% clf
		% plot(-maxP,max(Ssf),'+')
		% ylabel('fft power')
		% xlabel('-maxDepth')
		% grid on
		% title(sprintf('Individual %s Climb depth to fft power correlation',individual))
		% print([dir '/climb_depth2power.png'],'-dpng')

		% Correlation between swimming frequencies and maxP
		% clf
		% plot(-maxP,1./mean_f,'+')
		% ylabel('swim period (s)')
		% xlabel('-maxDepth')
		% grid on
		% title(sprintf('Individual %s Climb depth to swim period correlation',individual))
		% print([dir '/climb_depth2period.png'],'-dpng')

		% clf
		% plot3(t/86400,-maxP,max(Ssf),'.');
		% grid on
		% xlabel('t (j)')
		% ylabel('-maxP (m)')
		% zlabel('fft power')
		% title(sprintf('Individual %s Climb depth to fft power correlation, over time',individual))
		% print([dir '/climb_depthtime2power.png'],'-dpng')

		% Co-evolution of frequency power and max depth
		% clf
		% plot(t/86400,max(Ssf),'-k','linewidth',2)
		% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 60 10])
		% ylabel('max fft power')
		% xlabel('t (j)')
		% grid on
		% title(sprintf('Individual %s fft power evolution',individual))
		% print([dir '/climb_t2power.png'],'-dpng','-r100')

		% clf
		% plot(t/86400,std_f,'-k','linewidth',2)
		% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 60 10])
		% ylabel('std(f) (Hz)')
		% xlabel('t (j)')
		% grid on
		% title(sprintf('Individual %s std(f) evolution',individual))
		% print([dir '/climb_t2stdf.png'],'-dpng','-r100')

		% clf
		% plot(t/86400,-maxP,'-k')
		% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 60 10])
		% ylabel('-maxP (m)')
		% xlabel('t (j)')
		% grid on
		% title(sprintf('Individual %s dive depth evolution',individual))
		% print([dir '/climb_t2maxP.png'],'-dpng','-r100')

		% clf
		idx=find(-maxP>200);
		% plot(t(idx)/86400,max(Ssf(:,idx))./(-maxP(idx)),'.k')
		% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 60 10])
		% ylabel('fft power/maxP')
		% xlabel('t (j)')
		% grid on
		% title(sprintf('Individual %s normalised fft power evolution',individual))
		% print([dir '/climb_t2normpower.png'],'-dpng','-r100')

		is_good=zeros(size(t));
		is_good(idx)=ones(size(idx));
		climbfreqstats=[(t/86400)' is_good' -maxP' std_f' max(Ssf)'];
		% save([dir '/climbfreqstats.csv'],'-ascii','-double','climbfreqstats');

		% close all
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
		% clf
		% imagesc(t/86400,f(ifreq),Ssf)
		% xlabel('t (j)')
		% ylabel('f (Hz)')
		% title(sprintf('Individual %s Climb spectrogram',individual))
		% print([dir '/dive_spectro.png'],'-dpng')
		% Correlation between frequency power and max depth
		% clf
		% plot(-maxP,max(Ssf),'+')
		% ylabel('fft power')
		% xlabel('-maxDepth')
		% grid on
		% title(sprintf('Individual %s Climb depth to fft power correlation',individual))
		% print([dir '/dive_depth2power.png'],'-dpng')

		% Correlation between swimming frequencies and maxP
		% clf
		% plot(-maxP,1./mean_f,'+')
		% ylabel('swim period (s)')
		% xlabel('-maxDepth')
		% grid on
		% title(sprintf('Individual %s Climb depth to swim period correlation',individual))
		% print([dir '/dive_depth2period.png'],'-dpng')

		% clf
		% plot3(t/86400,-maxP,max(Ssf),'.');
		% grid on
		% xlabel('t (j)')
		% ylabel('-maxP (m)')
		% zlabel('fft power')
		% title(sprintf('Individual %s Climb depth to fft power correlation, over time',individual))
		% print([dir '/dive_depthtime2power.png'],'-dpng')

		% Co-evolution of frequency power and max depth
		% clf
		% plot(t/86400,max(Ssf),'-k','linewidth',2)
		% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 60 10])
		% ylabel('max fft power')
		% xlabel('t (j)')
		% grid on
		% title(sprintf('Individual %s fft power evolution',individual))
		% print([dir '/dive_t2power.png'],'-dpng','-r100')

		% clf
		% plot(t/86400,std_f,'-k','linewidth',2)
		% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 60 10])
		% ylabel('std(f) (Hz)')
		% xlabel('t (j)')
		% grid on
		% title(sprintf('Individual %s std(f) evolution',individual))
		% print([dir '/dive_t2stdf.png'],'-dpng','-r100')

		% clf
		% plot(t/86400,-maxP,'-k')
		% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 60 10])
		% ylabel('-maxP (m)')
		% xlabel('t (j)')
		% grid on
		% title(sprintf('Individual %s dive depth evolution',individual))
		% print([dir '/dive_t2maxP.png'],'-dpng','-r100')

		% clf
		idx=find(-maxP>200);
		% plot(t(idx)/86400,max(Ssf(:,idx))./(-maxP(idx)),'.k')
		% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 60 10])
		% ylabel('fft power/maxP')
		% xlabel('t (j)')
		% grid on
		% title(sprintf('Individual %s normalised fft power evolution',individual))
		% print([dir '/dive_t2normpower.png'],'-dpng','-r100')

		is_good=zeros(size(t));
		is_good(idx)=ones(size(idx));
		divefreqstats=[(t/86400)' is_good' -maxP' std_f' max(Ssf)'];
		% save([dir '/divefreqstats.csv'],'-ascii','-double','divefreqstats');

		% Desired aggregation:
		% Indiv, dive, min(tstart), max(tend), avg(drift), max(dP/dt),  sum(dive time), sum(climb time), max(max depth), 
		% avg(dive fft power), avg(climb fft power)
		% 
		diveid=unique(seg(2:end,1)); % Temp fix
		divestat=zeros(size(diveid,1),11);
		divestat(:,1) = ones(size(diveid,1),1)*individ;
		divestat(:,2) = diveid;
		for id=1:length(diveid)
			segi=seg(find(seg(:,1)==diveid(id)),:);
			divestat(id,3) = min(segi(:,7)/86400);
			divestat(id,4) = max(segi(:,9)/86400);
			if (divestat(id,4)-divestat(id,3))*86400 < 180.
				continue
			end
			idxdrift=find((drift(:,1)>=divestat(id,3)).*(drift(:,1)<=divestat(id,4)));
			if length(idxdrift)>0
				divestat(id,5) = mean(drift(idxdrift,2));
			end
			% idxpique=find((pique(:,1)>=divestat(id,3)).*(pique(:,1)<=divestat(id,4)));
			% if length(idxpique)>0
			% 	divestat(id,6) = max(pique(idxpique,2));
			% end
			idxdive=find(segi(:,5)==3);
			if length(idxdive)>0
				divestat(id,7) = sum((segi(idxdive,9)-segi(idxdive,7))/86400);
			end
			idxclimb=find(segi(:,5)==4);
			if length(idxclimb)>0
				divestat(id,8) = sum((segi(idxclimb,9)-segi(idxclimb,7))/86400);
			end
			idxdive=find((divefreqstats(:,1)>=divestat(id,3)).*(divefreqstats(:,1)<=divestat(id,4)));
			idxclimb=find((climbfreqstats(:,1)>=divestat(id,3)).*(climbfreqstats(:,1)<=divestat(id,4)));
			if length(idxdive)>0 && length(idxclimb)>0
				divestat(id,9) = max(max(divefreqstats(idxdive,3)),max(climbfreqstats(idxclimb,3)));
			end
			if length(idxdive)>0
				divestat(id,10) = mean(divefreqstats(idxdive,5));
			end
			if length(idxclimb)>0
				divestat(id,11) = mean(climbfreqstats(idxclimb,5));
			end
		end
		divestat(:,6) = interp1(pique(:,1),pique(:,2),(divestat(:,3)+divestat(:,4))/2,'previous');
		idxnan=find(isnan(divestat(:,6)));
		if length(idxnan)>0
			divestat(idxnan,6) = zeros(size(idxnan));
		end
		save([dir '/aggregate.csv'],'-ascii','-double','divestat');

		fid=fopen([dir '/individ.csv'],'w');
		fprintf(fid,'%d %s\n',individ,individual);
		fclose(fid);
		

		% disp('press enter to continue')
		% pause
	end
end
