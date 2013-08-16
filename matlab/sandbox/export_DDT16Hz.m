

matobj = matfile('DDT16hz_estimationvitesse')

batch=24*3600*16;
sz = size(matobj,'DDT');
i=1
cnt=1
txt=matobj.DDTtxt;
while  i <= sz(1)
    seekto=min(i+batch,sz(1));
    data=matobj.DDT(i:seekto,:);
    %filename=sprintf('DDT16Hz_%03d.txt',cnt);
    % save(filename,'-ascii','data')
    filename=sprintf('DDT16Hz_%03d.mat',cnt);
    save(filename,'data','txt');
    cnt = cnt + 1
    i = seekto+1
end


