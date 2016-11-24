
% load 'data/preload_mat_all.txt'
% P=preload_mat_all;
% clear preload_mat_all;
% load 'data/X'
% load 'data/Ap'
% load 'data/Bp'
% load 'data/Ap1'
% load 'data/Bp1'

n=size(Ap,1);
a=1
b=25000
figure(1);
plot(Ap(a:b,1),Ap(a:b,2),'r-',Ap(a:b,1),Ap(a:b,3),'g-',Ap(a:b,1),Ap(a:b,4),'b-')
grid on
figure(2);
plot(Bp(a:b,1),Bp(a:b,2),'r-',Bp(a:b,1),Bp(a:b,3),'g-',Bp(a:b,1),Bp(a:b,4),'b-')
% plot(Bp(a:b,1),Bp(a:b,2),'r-',Bp(a:b,1),Bp(a:b,3),'g-',Bp(a:b,1),Bp(a:b,4),'b-', Bp(a:b,1), P(a:b,9)/500, 'k-')
grid on
ab=[1 10e3;
    500e3 510e3;
    1000e3 1010e3];
for i=1:3
    a=ab(i,1);b=ab(i,2);
    figure(i)
    plot(Ap(a:b,1),Ap(a:b,2),'r-', 'linewidth',2,Ap(a:b,1),Ap(a:b,3),'g-', 'linewidth',2,Ap(a:b,1),Ap(a:b,4),'b-', 'linewidth',2,Ap(a:b,1),P(a:b,9)/200+5, 'k-','linewidth',3)
    grid on
    % print(sprintf('Ap%06d.png',a),'-dpng')
    % plot(Bp(a:b,1),Bp(a:b,2),'r-', 'linewidth',1,Bp(a:b,1),Bp(a:b,3),'g-', 'linewidth',1,Bp(a:b,1),Bp(a:b,4),'b-', 'linewidth',1,Bp(a:b,1),P(a:b,9)/2000+1, 'k-','linewidth',3)
    % grid on
    % % print(sprintf('Bp%06d.png',a),'-dpng')
    figure(3+i);
    plot(X(a:b,6),'r-', 'linewidth',3,P(a:b,9)/500-1, 'k-','linewidth',3)
    x = axis()
    x(3:4) = [-3 2];
    axis(x);
    grid on
    % print(sprintf('P%06d.png',a),'-dpng')
end

