function [d c]=kmccpdist(V,sB,subspace,A)
%function [d c]=kmccpdist(V,sB,subspace,A)
%calculates the distance of the points in
%sB to the plane V.  V(1,:) is the center of mass, 
%V(2:end,:) is a o.n. basis
%[tc tB]=kmccp(centers{k},B);
%d(:,k)=sum((tB-B).^2,2);
if nargin<3
    subspace=0;
end
[N,dim]=size(sB);
[cc,vdim]=size(V);
if nargin<4
    if subspace==1
        dt=sum(sB.^2,2);
        c=sB*V(1:cc,:)';
        dn=sum(c.^2,2);
        d=dt-dn;
    else
        tmcc=V(1,:);
        sB=sB-ones(N,1)*tmcc;
        dt=sum(sB.^2,2);
        c=sB*V(2:cc,:)';
        dn=sum(c.^2,2);
        d=dt-dn;
    end
else
    pB=(sB*V(1:cc,:)')*V(1:cc,:);
    dx=sB-pB;
    d=dx*A;
    d=sum(d.*dx,2);
end



