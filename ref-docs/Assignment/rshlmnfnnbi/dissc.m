function [mincent ttt]=dissc(dis,K,delta,opt) 
warning off all;
% dis: the N times N matrix representing the distances between N points and N subspaces
% K: the number of clusters
% opt: a N vector represeting the the noise around each points
% mincent: -- a structure of size length(delta)+2, each of the structure is a vector representing the classification arrording to the different delta values,
% and the last one is the one minimizing the ell_2 error
% ttt---a length(delta) times K+2 matrix, whose last row represent the eigenggap betwenn the squared K-th singular values and the squared K+1-th singular values.
dis=(dis.*dis').^0.5;
for j=1:size(dis,1)
    similarity2(j,:)=dis(j,:)/opt(j);
end
ttt=zeros(K+1,length(delta));
for i=1:length(delta)

similarity=exp(-(similarity2.^2)/delta(i));
 
similarity=(similarity+similarity')/2;

diag1=diag(max(sum(similarity')',0.00));

similarity=eye(size(similarity,1))+(diag1^-0.5*similarity*diag1^-0.5);

	[U,S]=svd(similarity,0);
%
U=U(:,1:K+1);
S=S(1:K+1,1:K+1);

eigenv=U(:,1:K)*(S(1:K,1:K)-eye(K))^0.5;


ttt(:,i)=[((diag(S)-1).^2)'];% reduce each singular values by one since eye(size(similarity,1)) was added before



[mincent{i},centers,temp{i}]= kmeans(eigenv,K,'start','sample','maxiter',800,'replicates',200*K,'EmptyAction','singleton');
end
N=size(ttt,1); 
for i=1:size(ttt,2)
ttt(N+1,i)=(ttt(N-1,i)^0.5-ttt(N,i)^0.5);

ttt(N+2,i)=(ttt(N-1,i)-ttt(N,i));% eigengap 
end
