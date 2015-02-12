function [Z e V m error]=pcp(X,N,subtractmean,density);

% function [Z e V m]=pcp(Z,N,subtractmean,density);

% function to subtract mean and project onto N pca vectors.  Z is

% samples x dim

% if subtractmean=0, data is not centered by pcp.  if subtractmean=1, pcp centers the

% data before doing principal components.

% density is a density on the data.  default is equal weight to every

% point. 

% output: Z=transformed data

% e=eigenvalues of covariance matrix

% V=eigenvectors of covariance matrix

% m=data mean 







if nargin==2

    subtractmean=1;

end

[samples dim]=size(X);





m=sum(X,1)/samples;

Z=X;

clear X;

if subtractmean==1;
%m=Z(1,:); 
for k=1:dim
	
    Z(:,k)=Z(:,k)-m(k);
	%Z(:,k)=Z(:,k)-Z(1,k);

end

end

%if dim<samples 
if 1==1
    if nargin<4

        q=Z'*Z;

    else

        W=Z;

        for j=1:samples

            W(j,:)=density(j)*W(j,:);

        end

        q=Z'*W;

    end

    q=full(q);

    [V,E]=eig(q);

    V=fliplr(V);

    e=diag(E);

    ee=sqrt(flipud(e));

    Z0=Z*V(:,1:N);
    
    error=norm(Z*V(:,N+1:end),'fro')/(samples^0.5);
    
    Z=Z0;

else

    q=Z*Z';

    q=full(q);

    [V,E]=eig(q);

    V=fliplr(V);

    V=V(:,1:N);

    %E=E(1:50,1:50)

    e=diag(E);

    e=sqrt(flipud(e));

    ee=e;

    e=e(1:N);

    E=diag(e);

    %V=real(V);

    %E=real(E);

    Zn=V*E;
    
    error=norm(V(:,N+1:end)*E,'fro')^2

    e=diag(E);

    R=V'*Z;

    R=diag(1./e)*R;

    V=R';

    Z=Zn;

end

e=ee;



