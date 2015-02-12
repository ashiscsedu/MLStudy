function [centers mincent error coeffs]=kmc2gnb(B,nc,cc,maxiter,subspace,quiet);
%[centers mincent error coeffs]=kmc2gnb(B,nc,cc,maxiter,subspace,quiet);

if nargin<6
    quiet=0;
end

if nargin<5
    subspace=0;
end

clear centers
clear error
clear errorc
[N,dim]=size(B);

rid=randpermute(N);

if nc==1
    if subspace==1
         [Z e V m]=pcp(B,cc,0);
         centers{1}=V(:,1:cc)';
         mincent=ones(N,1);
         error=sum(e(cc+1:end));
         coeffs=0;
         return
    end
end



 ste=4;
 scales=20;
for k=1:nc
    [uid centers{k}]=gnhbd(B,rid(k),cc,2*cc,ste,scales,subspace,1);
    %if subspace==0
    %    centers{k}=centers
end
nchanged=100;
d=zeros(N,nc);
it=0;
while nchanged>1
    it=it+1;
    if it>maxiter
        break
    end
    for k=1:nc
        if subspace==1
            d(:,k)=kmccpdist(centers{k},B,1);
        else            
            d(:,k)=kmccpdist(centers{k},B);
        end
    end
    [tmu mincent]=min(d,[],2);
    error=sum(tmu);
    if quiet==0
        error
    end
    if it>1
        nchanged=length(find(mincent-omincent));
        if quiet<1
            nchanged
        end
    end
    [mm idx]=mvi(mincent);
    omincent=mincent;
    %errorc(it)=error;

    for j=1:length(idx)
        if length(idx{j})>cc
            
            if subspace==1
                [Z e V m]=pcp(B(idx{j},:),cc,0);
                V=V(:,1:cc);
                centers{j}(1:cc,:)=V';
            else
                [Z e V m]=pcp(B(idx{j},:),cc,1);
                V=V(:,1:cc);
                centers{j}(1,:)=m; 
                centers{j}(2:cc+1,:)=V';
            end
        end
    end
    
end


