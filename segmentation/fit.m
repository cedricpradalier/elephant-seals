% dives=load('output/divesum.csv');
% dives=load('ml18_294dsens5.dives');
% seg=load('ml18_294dsens5.seg');
dives=load('ml17_301asens5.dives');
seg=load('ml17_301asens5.seg');


N=size(dives,1);
% idx=[1:N]';
idx=find((dives(:,2)>=1).*(dives(:,2)<=9));
t=dives(idx,1)/86400;
a=dives(idx,2);
P=dives(idx,3);

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

lt=linspace(min(t),max(t),60)';
la=linspace(0,10,100)';
[Gt,Ga]=meshgrid(lt,la);
gt=reshape(Gt,prod(size(Gt)),1);
ga=reshape(Ga,prod(size(Ga)),1);
% Am=[ones(size(gt)) gt ga ga.*ga gt.*gt ga.*gt];
% Am=[ones(size(gt)) gt ga ga.*ga];
Am=[ones(size(gt)) gt ga ];
gdP=Am*C2;
gdP=reshape(gdP,length(la),length(lt));
figure(1);
plot3(dives(:,1)/86400,dives(:,2),dives(:,3),'.b');
grid on
hold on; mesh(Gt,Ga,gdP); hold off
xlabel('t')
ylabel('az')
zlabel('dP/dt')

idx0=find(dives(:,2)<0.1);
t0=dives(idx0,1)/86400;
P0=dives(idx0,3);
A1=[ones(size(idx0)) t0 t0.*t0];
B=P0;
D1=A1\B;
err = A1*D1-B;
disp('A0: quadratic RMSE')
sqrt(sum(err.^2))
1 - mean(err.^2)/var(P)

A2=[ones(size(idx0)) t0 ];
B=P0;
D2=A2\B;
err = A2*D2-B;
disp('A0: linear RMSE')
sqrt(sum(err.^2))
1 - mean(err.^2)/var(P)

segd=seg(find(seg(:,4)==5),:);
ts=segd(:,6)/86400;
vs=(segd(:,9)-segd(:,7))./(segd(:,8)-segd(:,6));
figure(2)
% plot(dives(idx0,1)/86400,dives(idx0,3),'+b',lt,[ones(size(lt)) lt lt.*lt]*D1,'-+r')
subplot(2,1,1)
plot(dives(idx0,1)/86400,dives(idx0,3),'.b',lt,[ones(size(lt)) lt lt.*lt]*D1,'-+r',lt,[ones(size(lt)) lt ]*D2,'-+g');
grid on
axis([0 60 -3 0])
subplot(2,1,2)
plot(ts,vs,'*r');
grid on
axis([0 60 -1 0])

