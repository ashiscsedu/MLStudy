function C=cutva(bid,w,ty);
% function C=cutva(bid,sid,w,ty);
% ty is rc, nc, rcc, ncc, or cut
% let d=sum w
% vol(S)=sum_{j \in S} d
% cut is \sum_{i\in sid,j\in bid} w_{ij}
% rc= ratio cut= cut/|bid| + cut/|sid|  
% nc= normalized cut= cut/vol(bid)+cut/vol(sid)
% rcc= ratio cheeger cut= cut/min(|bid|,|sid|)
% ncc= normalized cheeger cut = cut/min(vol(bid),vol(sid))
% w needs to be symmetric

%if nargin==2
%    ty='rcc';
%end
%
%q=sum(w(:,sid),2);
%cut=sum(q(bid));
%
%if strcmp(ty,'cut')
%    C=cut;
%elseif strcmp(ty,'rc')
%    C=cut/length(sid)+cut/length(bid)
%elseif strcmp(ty,'nc')
%    vol=sum(w);
%    C=cut/sum(vol(sid))+cut/sum(vol(bid))
%elseif strcmp(ty,'rcc')
%    C=cut/min(length(sid),length(bid));
%elseif strcmp(ty,'ncc')
%    vol=sum(w);
%    C=cut/min(sum(vol(sid)),sum(vol(bid)));
%end
for i=1:size(bid)
w(bid{i},bid{i})=0;
end
C=max(max(w));
