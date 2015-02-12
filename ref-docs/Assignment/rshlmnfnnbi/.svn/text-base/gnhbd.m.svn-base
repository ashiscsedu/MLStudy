function [uid q r tdi noiseerror]=gnhbd(B,nnid,cc,startstep,step,scales,subspace,quiet,nB);

% find an appropriate neighborhood for any point by choosing the first  minimum
%input: B: the data matrix of size N*D, representing N points in D-dimensional space
% nnid: to find the nnid-th point's neighborhood
% cc: the dimension of the subspace
% starstep: the minimum size of the neighborhood:
% step: the size of neihghbrohoods increase by 'step' each time
% scales: the maximum number of different neighborhoods to try
% subspace: 0 for linear subspace, 1 for affine subspace
% quiet: set to 1 
% nB: the N-vector representing the squared norm of the N points
% output: uid: a vector, representing the indexes of the points that are in the fiited neighborhood
% q: the basis if the subspace
% r: the beta umber of the neighborhood
% tdi: the order of the N points, according to their distancec to the nnid-th point
% noiseerror: the averaged fitted l2 error in the neighborhood

% For more explanation of the code, see gnhbd2.m




if nargin<8

    quiet=1;

end

if nargin<9

    nB=0;

end



[N,dim]=size(B);



scales=min(scales,floor(N-startstep)/step);
if scales<2
dis=(sum((B-repmat(B(nnid,:),size(B,1),1)).^2,2))^0.5;


[td,tdi]=sort(dis);
    uid=tdi([1:startstep]);

     [Zu eu vu mu noiseerror]=pcp(B(uid,:),cc,1-subspace);


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




if length(nnid)==1

    d=zeros(N,1);

    if length(nB)==1

        nB=sum(B.^2,2);

    end

    dtb=sum(x.^2)+nB-2*B*x';



    [td tdi]=sort(dtb);

    nnid=tdi(1:startstep+step*scales);

end



clear q

or=inf;
 

for k=1:scales

    num=startstep+k*step;


    if num>N/2

        break

    end

    uid=(tdi(1:num));

    [Zu eu vu mu ]=pcp(B(uid,:),min(cc+5,dim),1-subspace);

    if subspace==0

        q(1,:)=mu';

        q(2:cc+1,:)=vu(:,1:cc)';

    else

        q(1:cc,:)=vu(:,1:cc)';

    end

    d=kmccpdist(q,B(uid,:),subspace);

    r(k)=sum(d)/(length(uid))/td(num);

    if r(k)>or

        break% when the first minimum is found

    end

    or=r(k);

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





