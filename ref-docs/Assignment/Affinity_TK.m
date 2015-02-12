function A=Affinity_TK(coord,D)


%Consider D as the kernel parameter (useful choices from 1:10)


Ntotal=size(coord,3);
frames=size(coord,2);


W=reshape(coord,2*frames,Ntotal); %measurement matrix
[U,S,V] = svd(W);

alpha=4;

 
X=V(:,1:D);
    
X0=diag(1./sqrt(sum(X.^2,2)))*X;
A=abs(X0*X0').^alpha;
A=A-diag(diag(A));

        
   
end






