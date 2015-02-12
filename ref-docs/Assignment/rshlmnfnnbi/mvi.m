function [b idx]=mvi(v);
%function [b idx]=mvi(v);
%finds unique values of v, stores in b;
%idx{k} is the indices of v corresponding to
%b(k)

[vv sv]=sort(v);
[b zz yy]=unique(vv,'first');
nN=length(b);

%b=vv(zz)
%vv=v(sv)
zz(nN+1)=length(v)+1;
idx{nN}=[];
for k=1:nN
    tt=[zz(k):zz(k+1)-1];
    idx{k}=sv(tt);
end