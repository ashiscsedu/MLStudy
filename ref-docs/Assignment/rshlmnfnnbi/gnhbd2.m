function [uid q r tdi noiseerror]=gnhbd2(B,nnid,cc,startstep,step,scales,subspace,quiet,nB);
% find an appropriate neighborhood for any point by choosing the first local minimum
%input: B: the data matrix of size N*D, representing N points in D-dimensional space
% nnid: to find the nnid-th point's neighborhood
% cc: the dimension of the subspace
% starstep: the minimum size of the neighborhood:
% step: the size of neihghbrohoods increase by 'step' each time
% scales: the maximum number of different neighborhoods to try
% subspace: 0 for linear subspace, 1 for affine subspace
% quiet: set to 1 if there is no output
% nB: the N-vector representing the squared norm of the N points
% output: uid: a vector, representing the indexes of the points that are in the fiited neighborhood
% q: the basis if the subspace
% r: the beta umber of the neighborhood
% tdi: the order of the N points, according to their distancec to the nnid-th point
% noiseerror: the averaged fitted l2 error in the neighborhood





if nargin<8

    quiet=1;% default: no output

end

if nargin<9

    nB=0;%we can leave the input of nB empty 

end



[N,dim]=size(B);



scales=min(scales,floor((N-startstep)/step));% the maxmimum steps is bounded by the number of points



if scales<2% if there is only one different neighborhoods to try, we just calculate the fitted subspace in the smallest neighborhood without calculating beta number

dis=(sum((B-repmat(B(nnid,:),size(B,1),1)).^2,2))^0.5;

[td,tdi]=sort(dis);
    uid=tdi([1:startstep]);

     [Zu eu vu mu noiseerror]=pcp(B(uid,:),cc,1-subspace);


    if subspace==0% for affine subspace, also calculate the center of mass

        q(1,:)=mu';

        q(2:cc+1,:)=vu(:,1:cc)';

    else

        q(1:cc,:)=vu(:,1:cc)';

    end

    r=0;

    tdi=0;

    return

end

x=B(nnid(1),:);%this is the point in the center of the neighborhood



    d=zeros(N,1);

    if length(nB)==1

        nB=sum(B.^2,2);

    end

    dtb=sum(x.^2)+nB-2*B*x';% squared distances of the N points to the center


    [td tdi]=sort(dtb);

    nnid=tdi(1:startstep+step*scales);





clear q
%tic;
test2=0;
test=0;
for k=1:scales

    num=startstep+k*step;

    uid=(nnid(1:num));

    [Zu eu vu mu]=pcp(B(uid,:),min(cc+5,dim),1-subspace);% fit a subspoace

    if subspace==0

        q(1,:)=mu';

        q(2:cc+1,:)=vu(:,1:cc)';% q is the description of the subspace

    else

        q(1:cc,:)=vu(:,1:cc)';

    end

    d=kmccpdist(q,B(uid,:),subspace);% distances of the points in the neiborhood to the fitted subspace

    r(k)=sum(d)/(length(uid)*td(num));%  the beta number for the k-th neighborhood
    
    if k==2
        test=((r(2)-r(1))>0);
    end
    
    if k>2
        if (r(k)>r(k-1) & r(k-1)<r(k-2))||(r(k)<r(k-1) & r(k-1)>r(k-2))
            test2=test2+1;
        end
    end
    if test2-test==1
        break;% when the first local minimum is found
    end

end

if k==scales% if we can not find the first local minimum, use the smallest neighborhood instead 
k=1;
end


num=startstep+(k-1)*step;% the size of the good neighborhood

uid=(tdi(1:num));% the label of the points in the good neighborhood

[Zu eu vu mu noiseerror]=pcp(B(uid,:),cc,1-subspace);% fit the good subspace




if subspace==0

    q(1,:)=mu';

    q(2:cc+1,:)=vu(:,1:cc)';

else

    q(1:cc,:)=vu(:,1:cc)';


end








