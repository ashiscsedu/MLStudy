function A=lbfsc(X,nc,cc,sigma)

planedimensions=ones(1,nc)*cc;
d=max(planedimensions);
K=length(planedimensions);
[N,D]=size(X);
step=2;% find the best neighborhood by increasing its size by 2 each time
scales=40;%check 40 different neighborhoods at most
ss=0;% affine subspace by default
nB=sum(X.^2,2);
rid=repmat(randpermute(floor(N)),1,1);
dis=zeros(N,N);
subspace=cell(1,N);
neighbor=cell(1,N);
neighbors=zeros(1,N);
delta=2*exp(0:1:6);% delta-- we multiple this parameter to the size of noise as the `delta' used in spectral clustering, could be a vector(default value: 2*exp(0:1:6))
for i=1:N 
	% fit a good subspace for each point

	[uid q q1 q2 noise(i)]=gnhbd(X,i,d,2*d,step,scales,ss,1,nB);

	subspace{i}=q;
	neighbor{i}=uid;
	neighbors(i)=length(uid);
	X_centered=X-repmat(subspace{i}(1,:),size(X,1),1);% reduce the center of the subspace from the data set
	dis(:,i)=(sum((X_centered-X_centered*subspace{i}(2:end,:)'*subspace{i}(2:end,:)).^2,2)).^0.5;%distance from the points to the subspaces

end	
	
A=exp(-dis.^2./sigma^2);