function [uid q r tdi noiseerror]=gnhbd2(B,nnid,cc,startstep,step,scales,subspace,quiet,nB);

%function [uid q r tdi]=gnhbd2(B,nnid,cc,startstep,step,scales,subspace,quiet,nB);





if nargin<8

    quiet=1;

end

if nargin<9

    nB=0;

end



[N,dim]=size(B);



scales=min(scales,floor((N-startstep)/step));



if scales<2
dis=(sum((B-repmat(B(nnid,:),size(B,1),1)).^2,2))^0.5;

%for i=1:N
%dis(i)=norm(B(i,:)-B(nnid,:));
%end
[td,tdi]=sort(dis);
    uid=tdi([1:startstep]);

     [Zu eu vu mu]=pcp(B(uid,:),min(cc+5,dim),1-subspace);

    if subspace==0

        q(1,:)=mu';

        q(2:cc+1,:)=vu(:,1:cc)';

    else

        q(1:cc,:)=vu(:,1:cc)';

    end

    r=0;

    tdi=0;

    return

end

x=B(nnid(1),:);
        for i=1:5
            sB=B-ones(N,dim)*diag(x);
            dt=sum(sB.^2,2);

            [td tdi]=sort(dt);

            x=mean(B(tdi(1:10),:),1);

        end


if length(nnid)==1

    d=zeros(N,1);

    if length(nB)==1

        nB=sum(B.^2,2);

    end

    dtb=sum(x.^2)+nB-2*B*x';

%    sB=B-ones(N,dim)*diag(x);

%    dtb=sum(sB.^2,2);

    %dtb=B*x';

    [td tdi]=sort(dtb);

    nnid=tdi(1:startstep+step*scales);

end



clear q
test2=0;
test=0;
for k=1:scales

    num=startstep+k*step;

    uid=(nnid(1:num));

    [Zu eu vu mu]=pcp(B(uid,:),min(cc+5,dim),1-subspace);

    if subspace==0

        q(1,:)=mu';

        q(2:cc+1,:)=vu(:,1:cc)';

    else

        q(1:cc,:)=vu(:,1:cc)';

    end

    d=kmccpdist(q,B(uid,:),subspace);

    %rs(k)=sum(d(2:end)./td(uid(2:end)))/(length(uid)-1);

    r(k)=sum(d)/(length(uid)*td(num));

 
    if k==2
        test=((r(2)-r(1))>0);
    end
    
    if k>2
        if (r(k)>r(k-1) & r(k-1)<r(k-2))||(r(k)<r(k-1) & r(k-1)>r(k-2))
            test2=test2+1;
        end
    end
    if test2-test==1
        break;
    end


end

if k==scales
k=1;
end



num=startstep+(k-1)*step;

uid=(tdi(1:num));

[Zu eu vu mu noiseerror]=pcp(B(uid,:),cc,1-subspace);


if subspace==0

    q(1,:)=mu';

    q(2:cc+1,:)=vu(:,1:cc)';

else

    q(1:cc,:)=vu(:,1:cc)';

end




if quiet==0

    r

    [Z e v m]=pcp(B,3,1);

    for k=1:scales

        y=zeros(N,1);

        num=startstep+k*step;

        uid=(tdi(1:num));

        y(uid)=1;

%        plotplot(Z',3,y);        

    end

    



end





